from pathlib import Path

from ligandcontacts import ConsensusLigandContacts, find_consensus_ligand_contacts

FIXTURES_DIR = Path(__file__).parent / "fixtures"

# CDK2_HUMAN(PDB 1GZ8)のうち、共結晶化リガンドMBP周辺(6A以内)の断片を抜粋した実データ
# ([fixtures/1gz8_mbp_fragment.pdb](fixtures/1gz8_mbp_fragment.pdb))。配列が不連続な断片なので
# 基準配列へのアラインメント位置自体は本来の生物学的な対応とは一致しない
# (このテストでは集計ロジック自体の検証に使い、対応するreference resnumの具体的な値は検証しない)。
FRAGMENT = FIXTURES_DIR / "1gz8_mbp_fragment.pdb"

# CDK20_HUMAN(UniProt Q8IZL9)の実配列。
CDK20_SEQUENCE = (
    "MDQYCILGRIGEGAHGIVFKAKHVETGEIVALKKVALRRLEDGFPNQALREIKALQEMED"
    "NQYVVQLKAVFPHGGGFVLAFEFMLSDLAEVVRHAQRPLAQAQVKSYLQMLLKGVAFCHA"
    "NNIVHRDLKPANLLISASGQLKIADFGLARVFSPDGSRLYTHQVATRWYRAPELLYGARQ"
    "YDQGVDLWSVGCIMGELLNGSPLFPGKNDIEQLCYVLRILGTPNPQVWPELTELPDYNKI"
    "SFKEQVPMPLEEVLPDVSPQALDLLGQFLLYPPHQRIAASKALLHQYFFTAPLPAHPSEL"
    "PIPQRLGGPAPKAHPGPPHIHDFHVDRPLEESLLNPELIRPFILEG"
)


def test_find_consensus_ligand_contacts_counts_repeated_structure():
    # 同じ構造を2回渡す(=同じ接触が2つの独立構造で再現されたことを模する)。
    # 最低要求カウントはmax(2, round(0.2 * 2)) = 2なので、両方に共通する残基は全て採用されるはず。
    result = find_consensus_ligand_contacts([FRAGMENT, FRAGMENT], CDK20_SEQUENCE)

    assert result.n_ligands == 2
    assert result.min_count == 2
    assert len(result.anchor_resnums) > 0
    assert set(result.anchor_resnums) == {r for r, c in result.contact_counts.items() if c >= 2}
    assert all(c == 2 for c in result.contact_counts.values())


def test_find_consensus_ligand_contacts_single_structure_below_min_count():
    # 1構造だけでは最低要求カウント(常に2以上)を満たせず、コンセンサス残基はゼロになる。
    result = find_consensus_ligand_contacts([FRAGMENT], CDK20_SEQUENCE)

    assert result.n_ligands == 1
    assert result.min_count == 2
    assert result.anchor_resnums == []
    assert len(result.contact_counts) > 0  # 生の接触自体は記録されている


def test_find_consensus_ligand_contacts_empty_input():
    result = find_consensus_ligand_contacts([], CDK20_SEQUENCE)

    assert result == ConsensusLigandContacts(anchor_resnums=[], n_ligands=0, contact_counts={}, min_count=0)
