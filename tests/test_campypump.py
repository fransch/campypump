"""Regression tests on operon regions from 8 PubMLST genomes with known cmeB status.

Expected identities are the Geneious Prime values from the Peru C. coli
spreadsheet (global alignment, 11168 and PA001 cmeB references).
"""
import glob
import os

import pytest

from campypump.core import Params, classify_genome, load_refs

DATA = os.path.join(os.path.dirname(__file__), "data")
EXPECTED = {  # genome prefix: (call, pid_RE_PA001, pid_WT_11168 or None)
    "id-48": ("wildtype-cmeB", None, 100.0),
    "id-150796": ("wildtype-cmeB", None, None),
    "id-151401": ("wildtype-cmeB", None, None),
    "id-105572": ("RE-cmeB", 100.0, None),
    "id-151430": ("RE-cmeB", None, None),
    "id-151856": ("wildtype-cmeB", 77.714, 89.978),   # Geneious values
    "id-152083": ("wildtype-cmeB", 78.172, 95.453),
    "id-152096": ("RE-cmeB", 97.438, 79.776),
}


@pytest.fixture(scope="module")
def refs():
    return load_refs()


@pytest.mark.parametrize("prefix", sorted(EXPECTED))
def test_calls_and_identities(prefix, refs):
    fa = glob.glob(os.path.join(DATA, f"{prefix}_*.fasta"))[0]
    row = classify_genome(fa, refs, Params())
    call, pid_re, pid_wt = EXPECTED[prefix]
    assert row["cmeB_type"] == call
    if pid_re is not None:
        assert row["cmeB_identity_RE_PA001_pct"] == pytest.approx(pid_re, abs=0.01)
    if pid_wt is not None:
        assert row["cmeB_identity_WT_11168_pct"] == pytest.approx(pid_wt, abs=0.01)
