from molstd import standardize_smiles


def test_standardize_strips_salt():
    # aspirin sodium salt -> aspirin (parent)
    result = standardize_smiles("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]")
    assert result == "CC(=O)Oc1ccccc1C(=O)O"


def test_standardize_invalid_smiles_returns_none():
    assert standardize_smiles("not a smiles") is None


def test_standardize_is_canonical_and_consistent():
    a = standardize_smiles("c1ccccc1O")
    b = standardize_smiles("Oc1ccccc1")
    assert a == b
