"""化合物のBemis-Murckoスキャフォールド計算。"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def compute_scaffold(smiles: str) -> str | None:
    """SMILESからBemis-Murckoスキャフォールド(canonical SMILES)を求める。パース失敗時はNone。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold_mol)
