import pytest

from kinasemotifs import find_kinase_motifs

# CDK20_HUMAN(UniProt Q8IZL9)の実配列。P-loop/触媒Lys/HRD/DFGを全て含む典型的なキナーゼドメイン。
CDK20_SEQUENCE = (
    "MDQYCILGRIGEGAHGIVFKAKHVETGEIVALKKVALRRLEDGFPNQALREIKALQEMED"
    "NQYVVQLKAVFPHGGGFVLAFEFMLSDLAEVVRHAQRPLAQAQVKSYLQMLLKGVAFCHA"
    "NNIVHRDLKPANLLISASGQLKIADFGLARVFSPDGSRLYTHQVATRWYRAPELLYGARQ"
    "YDQGVDLWSVGCIMGELLNGSPLFPGKNDIEQLCYVLRILGTPNPQVWPELTELPDYNKI"
    "SFKEQVPMPLEEVLPDVSPQALDLLGQFLLYPPHQRIAASKALLHQYFFTAPLPAHPSEL"
    "PIPQRLGGPAPKAHPGPPHIHDFHVDRPLEESLLNPELIRPFILEG"
)


def test_find_kinase_motifs_detects_all_motifs_in_cdk20():
    motifs = find_kinase_motifs(CDK20_SEQUENCE)

    assert motifs.p_loop == (11, 16)
    assert CDK20_SEQUENCE[motifs.p_loop[0] - 1 : motifs.p_loop[1]] == "GEGAHG"
    assert motifs.catalytic_lys == 33
    assert CDK20_SEQUENCE[motifs.catalytic_lys - 1] == "K"
    assert motifs.hrd == (125, 127)
    assert CDK20_SEQUENCE[motifs.hrd[0] - 1 : motifs.hrd[1]] == "HRD"
    assert motifs.dfg == (145, 147)
    assert CDK20_SEQUENCE[motifs.dfg[0] - 1 : motifs.dfg[1]] == "DFG"
    assert motifs.dfg_plus1 == 148


def test_anchor_resnums_combines_all_motifs():
    motifs = find_kinase_motifs(CDK20_SEQUENCE)

    anchors = motifs.anchor_resnums

    assert anchors == {11, 12, 13, 14, 15, 16, 33, 125, 126, 127, 145, 146, 147, 148}


def test_find_kinase_motifs_raises_when_sequence_has_no_kinase_motifs():
    with pytest.raises(ValueError):
        find_kinase_motifs("MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
