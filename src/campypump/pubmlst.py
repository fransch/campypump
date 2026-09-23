"""Retrieve isolate metadata (MLST ST / clonal complex + provenance) from the PubMLST API.

Uses BIGSdb_downloader (https://github.com/kjolley/BIGSdb_downloader) when a token directory is
given, so private records visible to your PubMLST account are included; otherwise anonymous
requests (public records only).

Route: GET /db/pubmlst_campylobacter_isolates/isolates/{id}
       -> provenance{...}, schemes[{description: "MLST", fields: {ST, clonal_complex}}], sequence_bin{...}
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

ISOLATE_URL = "https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates/isolates/{}"
SKIP_PROVENANCE = {"id", "sender", "curator", "datestamp"}
RENAME = {"species": "pubmlst_species"}   # 'species' is campypump's own --species label


def safe_name(key: str) -> str:
    """Stata-safe column name (letters, digits, underscore; <= 32 characters)."""
    k = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_")
    if k and k[0].isdigit():
        k = "v_" + k
    k = RENAME.get(k, k)
    return k[:32]


def fetch_record(iid: str, token_dir: str | None) -> tuple[dict | None, str]:
    """Return (json record or None, error message)."""
    url = ISOLATE_URL.format(iid)
    if token_dir:
        exe = shutil.which("bigsdb-downloader")
        if exe is None:
            return None, "bigsdb-downloader not installed (pip install bigsdb-downloader)"
        res = subprocess.run([exe, "--cron", "--key_name", "PubMLST", "--site", "PubMLST",
                              "--token_dir", token_dir, "--url", url],
                             capture_output=True, text=True)
        text = res.stdout
        if "{" not in text:
            return None, (res.stderr.strip().splitlines() or ["no response"])[-1]
        text = text[text.index("{"):]
    else:
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                text = resp.read().decode()
        except urllib.error.HTTPError as e:
            return None, f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            return None, str(e)
    try:
        rec = json.loads(text)
    except json.JSONDecodeError:
        return None, "response was not JSON"
    if "provenance" not in rec:
        return None, rec.get("message", "no provenance in response")
    return rec, ""


def parse_record(iid: str, rec: dict) -> dict:
    row: dict = {"pubmlst_id": iid}
    for k, v in rec.get("provenance", {}).items():
        if k in SKIP_PROVENANCE:
            continue
        row[safe_name(k)] = "; ".join(map(str, v)) if isinstance(v, list) else v
    if rec.get("aliases"):
        row["aliases"] = "; ".join(map(str, rec["aliases"]))
    sb = rec.get("sequence_bin") or {}
    row["assembly_contigs"] = sb.get("contig_count", "")
    row["assembly_length_bp"] = sb.get("total_length", "")
    for sch in rec.get("schemes", []):
        if str(sch.get("description", "")).strip().upper() == "MLST":
            for k, v in (sch.get("fields") or {}).items():
                row[safe_name(k)] = "; ".join(map(str, v)) if isinstance(v, list) else v
            row["mlst_loci_designated"] = sch.get("loci_designated_count", "")
    return row


def collect_metadata(ids: list[str], out_tsv: str, token_dir: str | None, delay: float = 0.5) -> None:
    """Fetch records for all ids (resumable via <out>.cache.jsonl) and write a TSV."""
    cache_path = out_tsv + ".cache.jsonl"
    cache: dict[str, dict] = {}
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    cache[str(r["pubmlst_id"])] = r
    failed = []
    todo = [i for i in ids if i not in cache]
    print(f"{len(ids)} ids: {len(ids) - len(todo)} already retrieved, {len(todo)} to fetch", file=sys.stderr)
    with open(cache_path, "a") as fh:
        for n, iid in enumerate(todo, 1):
            rec, err = fetch_record(iid, token_dir)
            if rec is None:
                failed.append((iid, err))
                print(f"[{n}/{len(todo)}] id {iid}: failed ({err})", file=sys.stderr)
            else:
                row = parse_record(iid, rec)
                cache[iid] = row
                fh.write(json.dumps(row) + "\n")
                fh.flush()
            if n % 100 == 0:
                print(f"[{n}/{len(todo)}] retrieved {n - len(failed)}, failed {len(failed)}  {time.strftime('%H:%M')}",
                      file=sys.stderr)
            time.sleep(delay)

    rows = [cache[i] for i in ids if i in cache]
    lead = ["pubmlst_id", "isolate", "pubmlst_species", "ST", "clonal_complex", "year", "source"]
    cols = [c for c in lead if any(c in r for r in rows)]
    cols += sorted({k for r in rows for k in r} - set(cols))
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols) + "\n")
    if failed:
        with open(out_tsv + ".failed.txt", "w") as fh:
            fh.writelines(f"{i}\t{e}\n" for i, e in failed)
    print(f"metadata for {len(rows)} of {len(ids)} isolates -> {out_tsv}"
          + (f"  ({len(failed)} failed, listed in {out_tsv}.failed.txt)" if failed else ""), file=sys.stderr)
