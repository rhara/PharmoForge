import pytest

from ligandprep import prepare_ligand_pdbqt


def test_prepare_ligand_pdbqt_writes_pdbqt_with_atoms(tmp_path):
    output_path = tmp_path / "ethanol.pdbqt"

    result = prepare_ligand_pdbqt("CCO", "ethanol", output_path)

    assert result == output_path
    text = output_path.read_text()
    assert text.startswith("REMARK") or "ATOM" in text
    assert "ATOM" in text


def test_prepare_ligand_pdbqt_rejects_invalid_smiles(tmp_path):
    with pytest.raises(ValueError):
        prepare_ligand_pdbqt("not a smiles", "bad", tmp_path / "bad.pdbqt")
