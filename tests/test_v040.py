"""Tests for 0.4.0: the four species reference genomes ship with the package and
campypump runs skani itself.

The point of 0.4.0 is that a user who installs campypump can get a species call without
finding four genomes and a second tool first, so most of what is worth testing is that
the genomes really are in the package and that the call degrades sensibly when skani is
not installed. skani itself is not exercised here: these tests must pass on a machine
that does not have it.
"""
import ast
import os
import subprocess

import pytest

from campypump import species as sp

EXPECTED = {
    "Cjejuni_NCTC11168": ("C. jejuni", 1_641_481),
    "Ccoli_RM4661": ("C. coli", 1_824_273),
    "Clari_RM2100": ("C. lari", 1_525_460),
    "Cupsaliensis_RM3940": ("C. upsaliensis", 1_628_950),
}


def test_four_species_reference_genomes_are_bundled():
    from importlib import resources
    d = resources.files("campypump").joinpath("data/species")
    names = {e.name for e in d.iterdir() if e.name.endswith(".fasta.gz")}
    assert names == {f"{stem}.fasta.gz" for stem in EXPECTED}


def test_bundled_genomes_unpack_complete_and_intact(tmp_path, monkeypatch):
    # a fresh cache, so this also covers the unpacking path rather than a warm directory
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.delenv(sp.SPECIES_REF_ENV, raising=False)
    files = sp.species_reference_files()
    assert len(files) == 4
    for path in files:
        stem = os.path.basename(path).rsplit(".", 1)[0]
        want_species, want_len = EXPECTED[stem]
        header = contigs = 0
        bases = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    contigs += 1
                    header += 1
                else:
                    bases += len(line.strip())
        assert contigs == 1, f"{stem} should be one complete contig"
        assert bases == want_len, f"{stem} is {bases} bp, expected {want_len}"
        assert sp.species_from_reference(path) == want_species


def test_unpacking_is_idempotent(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.delenv(sp.SPECIES_REF_ENV, raising=False)
    first = sp.species_reference_files()
    stamps = {f: os.stat(f).st_mtime_ns for f in first}
    second = sp.species_reference_files()
    assert first == second
    assert {f: os.stat(f).st_mtime_ns for f in second} == stamps, \
        "a warm cache must not be rewritten on every run"


def test_env_var_overrides_the_bundled_references(tmp_path, monkeypatch):
    own = tmp_path / "mine"
    own.mkdir()
    (own / "Ccoli_RM4661.fasta").write_text(">x\nACGT\n")
    monkeypatch.setenv(sp.SPECIES_REF_ENV, str(own))
    assert sp.species_reference_files() == [str(own / "Ccoli_RM4661.fasta")]


def test_skani_is_given_list_files_not_three_thousand_arguments(tmp_path, monkeypatch):
    """A Peru-sized run is ~3,700 genomes. Passing those paths on the command line would
    overflow the shell's argument limit, so query and reference paths go in list files."""
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.delenv(sp.SPECIES_REF_ENV, raising=False)
    genomes = []
    for i in range(3700):
        g = tmp_path / f"id-{i}.fasta"
        g.write_text(">c\nACGT\n")
        genomes.append(str(g))

    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        seen["ql"] = open(cmd[cmd.index("--ql") + 1]).read().splitlines()
        seen["rl"] = open(cmd[cmd.index("--rl") + 1]).read().splitlines()
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    sp.run_skani(genomes, str(tmp_path / "out.tsv"), threads=8)

    assert "--ql" in seen["cmd"] and "--rl" in seen["cmd"]
    assert len(seen["cmd"]) < 20, "paths must not be expanded onto the command line"
    assert len(seen["ql"]) == 3700
    assert len(seen["rl"]) == 4
    assert all(os.path.isabs(p) for p in seen["ql"] + seen["rl"])


def test_skani_failure_raises_with_the_command_shown(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.delenv(sp.SPECIES_REF_ENV, raising=False)
    g = tmp_path / "id-1.fasta"
    g.write_text(">c\nACGT\n")
    monkeypatch.setattr(sp.subprocess, "run",
                        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "boom"))
    with pytest.raises(RuntimeError, match="boom"):
        sp.run_skani([str(g)], str(tmp_path / "out.tsv"))


def test_missing_skani_leaves_typing_intact(tmp_path, monkeypatch, capsys):
    """Without skani the species columns are empty and the run still produces its calls:
    a missing optional dependency must not block the thing the tool is for."""
    from campypump import cli
    monkeypatch.setattr(sp, "skani_available", lambda: False)
    monkeypatch.setattr(cli, "sys", cli.sys)

    class A:
        no_ani = False
        ani = None
        output = str(tmp_path / "res.tsv")
        threads = 1

    rows = [{"sample": "id-1", "cmeB_type_final": "RE-cmeB"}]
    assert cli._add_species_call(A(), ["id-1.fasta"], rows) == 0
    assert rows[0]["cmeB_type_final"] == "RE-cmeB"
    assert "conda install -c bioconda skani" in capsys.readouterr().err


def test_no_ani_flag_skips_the_call(tmp_path, monkeypatch):
    from campypump import cli

    def explode(*_a, **_k):
        raise AssertionError("skani must not run when --no-ani is given")

    monkeypatch.setattr(sp, "run_skani", explode)

    class A:
        no_ani = True
        ani = None
        output = str(tmp_path / "res.tsv")
        threads = 1

    assert cli._add_species_call(A(), ["id-1.fasta"], [{"sample": "id-1"}]) == 0


def test_precomputed_table_is_used_instead_of_running_skani(tmp_path, monkeypatch):
    """The HPC path: a table made once for thousands of genomes, reused here."""
    from campypump import cli

    def explode(*_a, **_k):
        raise AssertionError("skani must not run when --ani is given")

    monkeypatch.setattr(sp, "run_skani", explode)
    table = tmp_path / "species_ani.tsv"
    table.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/r/Ccoli_RM4661.fasta\t/q/id-1.fasta\t98.1\t91.0\t90.4\tRM4661\tid-1\n"
        "/r/Cjejuni_NCTC11168.fasta\t/q/id-1.fasta\t86.2\t70.0\t69.1\tNCTC11168\tid-1\n")

    class A:
        no_ani = False
        ani = str(table)
        output = str(tmp_path / "res.tsv")
        threads = 1

    rows = [{"sample": "id-1", "species": "C. jejuni"}]
    flagged = cli._add_species_call(A(), ["/q/id-1.fasta"], rows)
    assert rows[0]["species_by_ANI"] == "C. coli"
    assert rows[0]["ANI_margin"] == "11.90"
    assert rows[0]["species_check_ANI"] == "CHECK"
    assert flagged == 1


def test_single_hit_leaves_the_margin_blank_without_failing():
    """A C. upsaliensis genome is below skani's reporting floor against the other three
    references, so it has no runner-up. That must not break the call."""
    import csv as _csv
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
                 "/r/Cupsaliensis_RM3940.fasta\t/q/id-9.fasta\t97.5\t88.0\t87.2\tRM3940\tid-9\n")
        path = fh.name
    try:
        rec = sp.read_ani(path)["id-9"]
    finally:
        os.unlink(path)
    assert rec["species_by_ANI"] == "C. upsaliensis"
    assert rec["ANI_margin"] == ""
    assert rec["ANI_call_quality"] == "ok"


def test_run_subcommand_exposes_both_ani_flags():
    """Guard the wiring: the flags must reach the parser, not merely exist in the code."""
    src = os.path.join(os.path.dirname(__file__), "..", "src", "campypump", "cli.py")
    tree = ast.parse(open(os.path.abspath(src)).read())
    flags = {a.value for node in ast.walk(tree)
             if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Attribute)
             and node.func.attr == "add_argument"
             for a in node.args if isinstance(a, ast.Constant) and isinstance(a.value, str)}
    assert "--ani" in flags and "--no-ani" in flags


# --------------------------------------------------------------------------- #
# 0.4.1: a directory argument is searched at every depth
# --------------------------------------------------------------------------- #

def test_nested_genome_directories_are_found(tmp_path):
    """genomes/<species>/id-N.fasta is the layout campypump's own fetch produces, so
    pointing run at the collection root must work."""
    from campypump.cli import _collect_genomes
    for sub, n in (("ccoli", 3), ("cjejuni", 2)):
        d = tmp_path / "genomes" / sub
        d.mkdir(parents=True)
        for i in range(n):
            (d / f"id-{sub}-{i}.fasta").write_text(">c\nACGT\n")
    found = _collect_genomes([str(tmp_path / "genomes")])
    assert len(found) == 5
    assert all(f.endswith(".fasta") for f in found)
    assert found == sorted(found)


def test_flat_directory_still_works(tmp_path):
    from campypump.cli import _collect_genomes
    for i in range(3):
        (tmp_path / f"id-{i}.fna").write_text(">c\nACGT\n")
    assert len(_collect_genomes([str(tmp_path)])) == 3


def test_empty_directory_says_where_it_looked(tmp_path, capsys):
    from campypump.cli import _collect_genomes
    (tmp_path / "notes.txt").write_text("nothing here")
    with pytest.raises(SystemExit) as e:
        _collect_genomes([str(tmp_path)])
    assert str(tmp_path) in str(e.value)
    assert ".fasta" in str(e.value)


def test_missing_path_is_reported_not_passed_to_blast(tmp_path):
    from campypump.cli import _collect_genomes
    with pytest.raises(SystemExit, match="No such file or directory"):
        _collect_genomes([str(tmp_path / "nope.fasta")])
