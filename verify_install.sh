#!/bin/bash
# ---------------------------------------------------------------------------
# campypump self-check. Run from the unpacked campypump folder:
#     bash verify_install.sh
# Runs the test suite, types the 8 bundled reference genomes, and compares the
# results with the values expected by the developers. Takes about 1 minute.
# ---------------------------------------------------------------------------
set -u
fail=0
out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT

echo "=== 1. Software ==="
python3 -c "import sys; print('python     ', sys.version.split()[0])" || fail=1
python3 -c "import Bio; print('biopython  ', Bio.__version__)" || fail=1
blastn -version 2>/dev/null | head -1 | sed 's/^/blast      /' || { echo "blastn NOT FOUND"; fail=1; }
campypump --version 2>/dev/null | sed 's/^/campypump  /' || { echo "campypump NOT INSTALLED"; fail=1; }

echo
echo "=== 2. Test suite ==="
python3 -m pytest -q tests || fail=1

echo
echo "=== 3. Typing the 8 bundled reference genomes ==="
campypump run tests/data/ -t 2 -o "$out/res.tsv" >/dev/null 2>"$out/log" || { cat "$out/log"; fail=1; }
cat "$out/log"

python3 - "$out/res.tsv" <<'PY'
import csv, sys

# expected values: cmeB type, % identity to PA001 RE-cmeB, % identity to NCTC 11168 cmeB.
# The three C. coli values are the laboratory's Geneious Prime results (tolerance 0.01).
EXPECTED = {
    "105572": ("RE-cmeB",       100.000,  79.598),
    "150796": ("wildtype-cmeB",  79.662,  97.983),
    "151401": ("wildtype-cmeB",  79.630,  98.591),
    "151430": ("RE-cmeB",        99.872,  79.534),
    "151856": ("wildtype-cmeB",  77.714,  89.978),   # Geneious
    "152083": ("wildtype-cmeB",  78.172,  95.457),   # Geneious: 95.453
    "152096": ("RE-cmeB",        97.438,  79.776),   # Geneious
    "48":     ("wildtype-cmeB",  79.598, 100.000),
}
rows = {r["pubmlst_id"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
bad = 0
print(f"{'id':8s} {'expected':14s} {'observed':14s} {'%RE':>8s} {'%WT':>8s}  result")
for iid, (call, pid_re, pid_wt) in sorted(EXPECTED.items()):
    r = rows.get(iid)
    if r is None:
        print(f"{iid:8s} {call:14s} {'MISSING':14s}"); bad += 1; continue
    got_call = r["cmeB_type_final"]
    got_re, got_wt = float(r["cmeB_identity_RE_PA001_pct"]), float(r["cmeB_identity_WT_11168_pct"])
    ok = got_call == call and abs(got_re - pid_re) < 0.01 and abs(got_wt - pid_wt) < 0.01
    bad += not ok
    print(f"{iid:8s} {call:14s} {got_call:14s} {got_re:8.3f} {got_wt:8.3f}  {'OK' if ok else 'MISMATCH'}")
sys.exit(1 if bad else 0)
PY
[ $? -ne 0 ] && fail=1

echo
if [ $fail -eq 0 ]; then
    echo "ALL CHECKS PASSED"
else
    echo "SOME CHECKS FAILED (see above)"
fi
exit $fail
