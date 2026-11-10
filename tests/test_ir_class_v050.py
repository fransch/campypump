"""0.5.0: the CmeR operator is classified by structure, not by identity to one strain.

CmeR binds as a dimer, one subunit per half-site, so what matters is whether the two
half-sites still pair and whether the spacer still holds them the right distance apart.
NCTC 11168 and 81-176 differ by one base in this window and are both intact operators;
classifying against either one alone makes the other look altered.
"""
import pytest

from campypump.core import classify_ir

CASES = [
    # sequence,                 half, spacer, length, class
    ("TGTAATAAATATTACA",          7,  2, 16, "intact"),                       # 81-176
    ("TGTAATAAAAATTACA",          6,  4, 16, "intact"),                       # NCTC 11168
    ("TGTAATAATATTACA",           7,  1, 15, "spacer_changed"),               # Lin 2005 lesion
    ("TGTAATAATAAATATTACA",       7,  5, 19, "spacer_changed"),
    ("TGTAATAAATATTACTATTACA",    7,  8, 22, "spacer_changed"),
    ("TGTAATAAATATTATA",          0, "", 16, "half_site_changed"),
    ("TGCAATAAATATTACA",          0, "", 16, "half_site_changed"),
    ("TGCAATAATAAATATTACA",       0, "", 19, "half_site_and_spacer_changed"),
]


@pytest.mark.parametrize("seq,half,spacer,length,cls", CASES)
def test_classify_ir(seq, half, spacer, length, cls):
    assert classify_ir(seq) == (half, spacer, length, cls)


def test_both_reference_strains_are_intact():
    """The whole point: neither reference is privileged."""
    assert classify_ir("TGTAATAAATATTACA")[-1] == classify_ir("TGTAATAAAAATTACA")[-1] == "intact"


def test_spacer_blank_when_half_sites_broken():
    """A spacer is undefined once the half-sites no longer pair."""
    assert classify_ir("TGTAATAAATATTATA")[1] == ""


@pytest.mark.parametrize("seq", ["", "TGTAAT", "ACGT"])
def test_short_windows_are_unclassified(seq):
    assert classify_ir(seq)[-1] == "unclassified"
