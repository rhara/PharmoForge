from molstd import calc_mol_weight


def test_calc_mol_weight_aspirin():
    # aspirin: C9H8O4, MW ~180.16
    mw = calc_mol_weight("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert mw is not None
    assert round(mw, 2) == 180.16


def test_calc_mol_weight_invalid_smiles_returns_none():
    assert calc_mol_weight("not a smiles") is None
