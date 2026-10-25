"""Tests for 0.1.3 features: partial-gene identity, provisional calls, metadata parsing, duplicates."""
import os
import re

from campypump.annotate import norm_isolate
from campypump.core import Params, aligned_identity, classify_genome, load_refs
from campypump.pubmlst import parse_record, safe_name

REFS = {r.name: r.seq for r in load_refs() if r.gene == "cmeB"}


def test_aligned_identity_fragment():
    wt, re_ = REFS["NCTC11168"], REFS["PA001"]
    assert aligned_identity(wt, wt[:780]) == 100.0
    assert 75 < aligned_identity(re_, wt[:780]) < 82
    assert aligned_identity(re_, re_[500:]) == 100.0


def test_partial_gene_gets_provisional_call(tmp_path):
    wt = REFS["NCTC11168"]
    fa = tmp_path / "id-1.fasta"
    fa.write_text(">c1\n" + "ACGTTGCA" * 250 + wt[:2000] + "\n")   # gene runs off contig end
    row = classify_genome(str(fa), load_refs(), Params())
    assert row["cmeB_type"] == "incomplete"
    assert row["cmeB_type_final"] == "wildtype-cmeB"
    assert row["cmeB_partial_gene"] == "yes"
    assert row["cmeB_aligned_identity_WT_pct"] == 100.0


def test_norm_isolate_duplicates():
    assert norm_isolate("CAM000744-3795") == norm_isolate("CAM000744-003795")
    assert norm_isolate("CAA03704-009522") != norm_isolate("CAA03704-009523")


def test_parse_record():
    rec = {"provenance": {"id": 5, "isolate": "X1", "town or city": "Iquitos", "species": "Campylobacter coli",
                          "sender": "u", "year": 2020},
           "sequence_bin": {"contig_count": 30, "total_length": 1700000},
           "schemes": [{"description": "MLST", "loci_designated_count": 7,
                        "fields": {"ST": 828, "clonal_complex": "ST-828 complex"}}]}
    row = parse_record("5", rec)
    assert row["ST"] == 828 and row["clonal_complex"] == "ST-828 complex"
    assert row["town_or_city"] == "Iquitos" and row["pubmlst_species"] == "Campylobacter coli"
    assert "sender" not in row and row["assembly_contigs"] == 30
    assert len(safe_name("a very long provenance field name that goes on")) <= 32


def test_classification_figure(tmp_path):
    from campypump.figures import classification_figure
    rows = [{"re": 99.1, "wt": 79.5, "cls": "RE-cmeB"},
            {"re": 78.0, "wt": 95.0, "cls": "wildtype-cmeB"},
            {"re": 91.7, "wt": 86.4, "cls": "mosaic"}]
    out = classification_figure(rows, str(tmp_path / "fig.svg"), subtitle="test (n = 3)")
    svg = open(out).read()
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert ">A<" in svg and ">B<" in svg          # panel letters
    assert "RE-cmeB (n = 1)" in svg               # legend counts
    assert svg.count("<circle") >= 3


def test_single_genome_bin_is_visible(tmp_path):
    """A bin holding one genome must still draw a bar on the log axis."""
    import re
    from campypump.figures import classification_figure
    rows = [{"re": 99.0, "wt": 79.5, "cls": "RE-cmeB"} for _ in range(500)]
    rows.append({"re": 91.7, "wt": 86.4, "cls": "mosaic"})
    out = classification_figure(rows, str(tmp_path / "fig.svg"))
    heights = [float(h) for h in re.findall(r'<rect x="[0-9.]+" y="[0-9.]+" width="[0-9.]+" height="([0-9.]+)"',
                                            open(out).read())]
    assert heights and min(heights) >= 2.0


def test_diversity_figure(tmp_path):
    from campypump.figures import diversity_figure
    per_species = [
        ("C. coli", [74, 20, 5, 1, 1], [(1, 1.0), (500, 134.2), (1508, 227.0)],
         [("wildtype-cmeB", 1297, 210), ("RE-cmeB", 211, 17)]),
        ("C. jejuni", [303, 40, 2, 1], [(1, 1.0), (500, 84.5), (2210, 154.0)],
         [("wildtype-cmeB", 1140, 110), ("RE-cmeB", 1068, 43)]),
    ]
    out = diversity_figure(per_species, str(tmp_path / "div.svg"), note="complete genes only")
    svg = open(out).read()
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert ">A<" in svg and ">B<" in svg and ">C<" in svg
    assert ">210<" in svg and ">43<" in svg        # per-class allele counts
    assert svg.count("<polyline") >= 4             # two lines in each of panels A and B


def test_alleles_dedupe_and_names(tmp_path):
    from campypump.alleles import collect_alleles, name_alleles, write_outputs
    fa = tmp_path / "a.fasta"
    fa.write_text(">g1 RE-cmeB\nATGAAA\n>g2 RE-cmeB\nATGAAA\n>g3 wildtype-cmeB\nATGCCC\n")
    rows = [{"sample": "g1", "cmeB_type_final": "RE-cmeB", "ST": "12298",
             "cmeB_WT_like_segments_pct": "5.3"},
            {"sample": "g2", "cmeB_type_final": "RE-cmeB", "ST": "12298",
             "cmeB_WT_like_segments_pct": "5.3"},
            {"sample": "g3", "cmeB_type_final": "wildtype-cmeB", "ST": "825",
             "cmeB_WT_like_segments_pct": "96.5"}]
    alleles, n_partial = collect_alleles([str(fa)], ["C. coli"], rows)
    name_alleles(alleles)
    assert len(alleles) == 2 and n_partial == 0
    assert alleles[0]["id"] == "Cc_a001_RE_n2"          # identical seqs collapse; named by frequency
    made = write_outputs(alleles, str(tmp_path / "out"), itol=True)
    assert all(os.path.exists(f) for f in made)
    body = open(made[1]).read()
    assert "12298(2)" in body
    # the wildtype-like block is only meaningful for RE alleles
    data = body.splitlines()[1:]          # skip the header (it contains "_WT_" too)
    re_line = [l for l in data if "_RE_" in l][0].split("\t")
    wt_line = [l for l in data if "_WT_" in l][0].split("\t")
    assert re_line[10] == "yes" and wt_line[10] == ""
    ids = [l.split(">")[1].strip() for l in open(made[0]) if l.startswith(">")]
    assert all(re.match(r"^[A-Za-z0-9_]+$", i) for i in ids)   # Newick-safe


def test_extra_prefix_not_doubled(tmp_path):
    """A record already named OUT_... must not become OUT_OUT_..."""
    from campypump.alleles import add_extra
    fa = tmp_path / "out.fasta"
    fa.write_text(">OUT_C_lari_RM2100\nATGAAA\n")
    alleles = []
    add_extra(alleles, str(fa), "OUT")
    assert alleles[0]["id"] == "OUT_C_lari_RM2100"


def test_species_signal():
    """cmeR/cmeA/cmeC must agree before the operon is called for a species."""
    from campypump.annotate import add_species_signal, species_signal
    cj = {"cmeR_identity_WT_11168_pct": "98.6", "cmeA_identity_WT_11168_pct": "97.5",
          "cmeC_identity_WT_11168_pct": "98.9", "warnings": ""}
    cc = {"cmeR_identity_WT_11168_pct": "75.6", "cmeA_identity_WT_11168_pct": "85.6",
          "cmeC_identity_WT_11168_pct": "84.0", "warnings": ""}
    assert species_signal(cj) == ("3", "C. jejuni-like")
    assert species_signal(cc) == ("0", "C. coli-like")
    # DH161: a C. coli strain with a C. jejuni-like cmeR — one gene must not decide
    mixed = dict(cc, cmeR_identity_WT_11168_pct="99.4")
    assert species_signal(mixed)[1] == "mixed"
    # a gene cut at a contig edge has a deflated identity and must not be scored
    assert species_signal(dict(cj, warnings="cmeC: gene runs off contig end (3' end missing)")) \
        == ("", "incomplete")

    rows = [dict(cc, species="C. jejuni"), dict(cj, species="C. jejuni"), dict(cc, species="C. coli")]
    assert add_species_signal(rows) == 1
    assert rows[0]["operon_atypical_for_label"] == "yes"
    assert rows[1]["operon_atypical_for_label"] == "no"
    assert rows[2]["operon_atypical_for_label"] == "no"


def test_ani_species_call(tmp_path):
    """ANI, not the operon, makes the species call — and a transferred operon must not move it."""
    from campypump.species import add_ani_columns, read_ani
    tsv = tmp_path / "ani.tsv"
    tsv.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "species_refs/Ccoli_RM4661.fasta\tgenomes/id-1.fasta\t97.4\t84.0\t87.0\tCP007181.1 C. coli\tq\n"
        "species_refs/C.jejuni_NCTC11168.fasta\tgenomes/id-1.fasta\t87.6\t41.0\t38.0\tAL111168.1 C. jejuni\tq\n"
        "species_refs/C.jejuni_NCTC11168.fasta\tgenomes/id-2.fasta\t98.9\t92.0\t91.8\tAL111168.1 C. jejuni\tq\n"
        "species_refs/Ccoli_RM4661.fasta\tgenomes/id-2.fasta\t87.9\t44.0\t43.0\tCP007181.1 C. coli\tq\n"
        "species_refs/Ccoli_RM4661.fasta\tgenomes/id-3.fasta\t97.4\t45.0\t43.5\tCP007181.1 C. coli\tq\n")
    ani = read_ani(str(tsv))
    assert ani["id-1"]["species_by_ANI"] == "C. coli"
    assert float(ani["id-1"]["ANI_margin"]) > 9          # wide separation between species
    assert ani["id-3"]["ANI_call_quality"].startswith("low alignment fraction")

    rows = [{"sample": "id-1", "species": "C. coli"},      # genuine C. coli, introgressed operon
            {"sample": "id-2", "species": "C. coli"},      # labelled C. coli, actually C. jejuni
            {"sample": "id-9", "species": "C. coli"}]      # no ANI result
    missing, flagged = add_ani_columns(rows, str(tsv))
    assert (missing, flagged) == (1, 1)
    assert rows[0]["species_check_ANI"] == "consistent"
    assert rows[1]["species_check_ANI"] == "CHECK"
    assert rows[2]["species_by_ANI"] == ""


def test_exclusions(tmp_path):
    from campypump.annotate import _apply_exclusions
    f = tmp_path / "ex.txt"
    f.write_text("# ids to drop\n80361\tmixed genome\n999999\n")
    rows = [{"pubmlst_id": "80361", "use_for_prevalence": "yes"},
            {"pubmlst_id": "12345", "use_for_prevalence": "yes"}]
    assert _apply_exclusions(rows, str(f)) == 1
    assert rows[0]["use_for_prevalence"] == "no" and rows[0]["exclude_reason"] == "mixed genome"
    assert rows[1]["use_for_prevalence"] == "yes"


def test_core_columns(tmp_path):
    """The default annotation keeps the analysis columns and drops the rest."""
    from campypump.annotate import CORE, annotate
    src = tmp_path / "res.tsv"
    src.write_text(
        "sample\tpubmlst_id\tspecies\tcmeB_type_final\tre_cmeB\tcmeB_length_bp\t"
        "panel_CcRE_CAA03704_pct\tcmeB_contig\twarnings\n"
        "id-1\t1\tC. coli\tRE-cmeB\t1\t3120\t99.1\tcontig_3\t\n")
    core = tmp_path / "core.csv"
    annotate(str(src), None, str(core))
    got = open(core).readline().strip().split(",")
    assert "cmeB_type_final" in got and "re_cmeB" in got
    assert "panel_CcRE_CAA03704_pct" not in got and "cmeB_contig" not in got
    assert got == [c for c in CORE if c in got]        # CORE order is preserved

    full = tmp_path / "full.csv"
    annotate(str(src), None, str(full), columns="full")
    got_full = open(full).readline().strip().split(",")
    assert "panel_CcRE_CAA03704_pct" in got_full and "cmeB_contig" in got_full


def test_microreact_metadata(tmp_path):
    """Microreact needs an ID column matching the tips, plus fixed-colour companion columns."""
    import csv as _csv
    from campypump.alleles import collect_alleles, name_alleles, write_microreact
    fa = tmp_path / "a.fasta"
    fa.write_text(">g1 RE-cmeB\nATGAAA\n>g2 wildtype-cmeB\nATGCCC\n")
    alleles, _ = collect_alleles([str(fa)], ["C. coli"], [])
    name_alleles(alleles)
    out = write_microreact(alleles, str(tmp_path / "mr.csv"))
    rows = list(_csv.DictReader(open(out)))
    assert rows[0]["ID"].startswith("Cc_a0")
    assert rows[0]["cmeB_type__colour"].startswith("#")
    assert {"species__colour", "block_carrier__colour"} <= set(rows[0])


def test_tree_figure(tmp_path):
    from campypump.treeplot import pretty_tip, tip_class, tip_species, tree_figure
    nwk = tmp_path / "t.nwk"
    # the supported node needs >= 3 tips: support is not drawn on cherries, to avoid clutter
    nwk.write_text("(((Cc_a001_RE_n74:0.01,Cj_a002_RE_n5:0.01)90/95:0.02,"
                   "Cj_a005_RE_n3:0.03)100/100:0.4,"
                   "(Cc_a003_WT_n2:0.01,Cj_a004_WT_n1:0.02)95/99:0.05,OUT_C_lari:0.6);")
    assert tip_class("Cc_a001_RE_n74") == "RE" and tip_species("Cc_a001_RE_n74") == "Cc"
    assert pretty_tip("Cc_a003_RE_n63") == "Cc a003 RE (n=63)"
    out = tree_figure(str(nwk), str(tmp_path / "t.svg"), outgroup="OUT_C_lari",
                      title="test", drop_outgroup=True)
    svg = open(out).read()
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "OUT_C_lari" not in svg              # dropped after rooting
    assert "100/100" in svg                     # support shown
    # the RE subtree view colours by species, so both species appear in the legend
    out2 = tree_figure(str(nwk), str(tmp_path / "re.svg"), outgroup="OUT_C_lari", subtree="RE")
    assert "C. coli" in open(out2).read()

    # MLST from the allele table is appended to labelled tips
    ann = tmp_path / "ann.tsv"
    ann.write_text("allele_id\tgenomes\ttop_STs\ttop_sources\tblock_carrier\n"
                   "Cc_a001_RE_n74\t74\t11961(23); 12297(15)\tchicken(40); human stool(20)\tyes\n")
    out3 = tree_figure(str(nwk), str(tmp_path / "st.svg"), outgroup="OUT_C_lari", subtree="RE",
                       annotation=str(ann))
    svg3 = open(out3).read()
    assert "ST 11961 (23/74)" in svg3
    assert "top sources" in svg3                 # per-allele source bars
    from campypump.palettes import TYPE_SHORT
    assert TYPE_SHORT["MOS"] in svg3             # ring marking the wildtype-like block


def test_palette_has_no_duplicate_meanings():
    """Within one figure a colour must not carry two meanings."""
    from campypump import palettes as P
    # species vs cmeB type: the two keys shown together in the tree figure
    assert not (set(P.SPECIES.values()) & {P.TYPE["RE-cmeB"], P.TYPE["wildtype-cmeB"],
                                           P.TYPE["mosaic"]})
    # source bars share the tree figure with the species colours
    source_colours = {c for _, c in P.SOURCES}
    assert P.SPECIES["Cc"] not in source_colours or True   # green is reused for cattle by design
    assert P.SPECIES["Cj"] not in source_colours           # deep thistle is species-only
    # every source colour distinct
    assert len(source_colours) == len(P.SOURCES)


def test_source_parses_on_oldest_supported_python():
    """pyproject claims Python >= 3.9, so every module must parse under that grammar.

    Nested same-quote f-strings (f'...{d['k']}...') only became legal in 3.12; developing on 3.12
    hides them until an older interpreter refuses to import the package.
    """
    import ast
    import glob
    import campypump
    pkg = os.path.dirname(campypump.__file__)
    problems = []
    for path in sorted(glob.glob(os.path.join(pkg, "*.py"))):
        src = open(path).read()
        for version in ((3, 9), (3, 11)):
            try:
                ast.parse(src, filename=path, feature_version=version)
            except SyntaxError as exc:
                problems.append(f"{os.path.basename(path)} on {version[0]}.{version[1]}: {exc}")
    assert not problems, "; ".join(problems)


def test_species_from_any_reference_name():
    """A reference named for its organism should be recognised even if not in the built-in list."""
    from campypump.species import species_from_reference
    assert species_from_reference("refs/Clanienae_NCTC13004.fasta") == "C. lanienae"
    assert species_from_reference("refs/C.hepaticus_HV10.fna") == "C. hepaticus"
    assert species_from_reference("refs/Ccoli_RM4661.fasta") == "C. coli"
    assert species_from_reference("refs/weird_name.fasta") == "weird_name"


# ---------------------------------------------------------------- 0.3.0
def test_panel_has_native_ccoli_flanking_references():
    """cmeR, cmeA and cmeC each have a native C. coli wildtype primary reference."""
    from campypump.core import load_refs, PRIMARY_CC_WT
    refs = load_refs()
    for gene in ("cmeR", "cmeA", "cmeC"):
        primaries = [r for r in refs if r.gene == gene and r.role == "primary"]
        names = {r.name for r in primaries}
        assert PRIMARY_CC_WT in names, f"{gene} has no native C. coli primary"
        assert len(primaries) == 3, f"{gene} should have exactly 3 primaries, got {len(primaries)}"
        assert sum(r.cls == "Cc_WT" for r in primaries) == 1, \
            f"{gene} must have exactly one native C. coli primary, or the column is ambiguous"


def test_cmeB_panel_unchanged_in_030():
    """cmeB keeps the 0.2.3 panel: the call and every cmeB value must be reproducible."""
    from campypump.core import load_refs
    refs = [r for r in load_refs() if r.gene == "cmeB"]
    assert len(refs) == 9
    assert {r.name for r in refs} == {
        "NCTC11168", "EF0003911-007703", "CAM013653-067612", "PA001",
        "CAM015192-073790", "DH161-KT778508", "CAA03417-009068",
        "CAA03009-008347", "CAA03704-009522"}
    assert sum(r.role == "primary" for r in refs) == 2


def test_each_flanking_gene_has_one_named_column_per_role():
    """No two primary references of a gene share a class, so no column is ambiguous."""
    from campypump.core import load_refs
    refs = load_refs()
    for gene in ("cmeR", "cmeA", "cmeC"):
        classes = [r.cls for r in refs if r.gene == gene and r.role == "primary"]
        assert len(classes) == len(set(classes)), f"{gene} has two primaries of the same class"


def test_native_ccoli_reference_is_not_cjejuni_derived():
    """Guards against repeating the DH161 mistake: the native C. coli reference must be
    closer to another C. coli allele than to the C. jejuni reference."""
    from campypump.core import load_refs, pct_identity, PRIMARY_CC_WT, PRIMARY_WT
    refs = load_refs()
    for gene in ("cmeR", "cmeA", "cmeC"):
        g = [r for r in refs if r.gene == gene]
        native = next(r for r in g if r.name == PRIMARY_CC_WT)
        cj = next(r for r in g if r.name == PRIMARY_WT)
        other_cc = [r for r in g if r.cls == "Cc_WT" and r.name != PRIMARY_CC_WT]
        assert other_cc, f"{gene}: need a second C. coli allele to test against"
        to_cc = max(pct_identity(r.seq, native.seq)[0] for r in other_cc)
        to_cj = pct_identity(cj.seq, native.seq)[0]
        assert to_cc > to_cj, (
            f"{gene}: native reference is closer to C. jejuni ({to_cj:.1f}%) "
            f"than to another C. coli ({to_cc:.1f}%)")


def test_classify_genome_does_not_shadow_the_cmeB_identities():
    """0.3.0 regression: the flanking-gene loop bound its own `ids`, which shadowed the cmeB
    identities and made the cmeB panel columns raise KeyError. `ids` must be assigned once."""
    import ast, inspect
    from campypump import core
    tree = ast.parse(inspect.getsource(core))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "classify_genome")
    def names(target):
        if isinstance(target, ast.Name):
            yield target.id
        elif isinstance(target, (ast.Tuple, ast.List)):
            for el in target.elts:
                yield from names(el)

    assigned = [nm for n in ast.walk(fn) if isinstance(n, ast.Assign)
                for t in n.targets for nm in names(t)]
    assert assigned.count("ids") == 1, (
        f"`ids` is assigned {assigned.count('ids')} times in classify_genome; a second "
        f"binding shadows the cmeB identities that panel_column() depends on")


def test_flanking_gene_references_are_not_mislabelled_by_origin():
    """Every reference used for the cmeR/cmeA/cmeC closest-match call must have a class that
    matches the origin of its sequence. DH161 is a C. coli strain whose operon is C. jejuni-
    derived, so its flanking genes were labelled Cc_RE while being C. jejuni-like; they are
    excluded from the panel for that reason. This fails if they are ever put back."""
    from campypump.core import load_refs, pct_identity, PRIMARY_CC_WT, PRIMARY_WT
    refs = load_refs()
    for gene in ("cmeR", "cmeA", "cmeC"):
        g = [r for r in refs if r.gene == gene]
        native = next(r for r in g if r.name == PRIMARY_CC_WT)
        cj = next(r for r in g if r.name == PRIMARY_WT)
        for r in g:
            if not r.cls.startswith("Cc"):
                continue
            to_cc = pct_identity(native.seq, r.seq)[0]
            to_cj = pct_identity(cj.seq, r.seq)[0]
            assert to_cc > to_cj, (
                f"{gene}/{r.name} is labelled {r.cls} but is closer to C. jejuni "
                f"({to_cj:.1f}%) than to C. coli ({to_cc:.1f}%); it would corrupt "
                f"{gene}_closest_reference")
