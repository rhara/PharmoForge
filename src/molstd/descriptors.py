"""RDKitを用いた化合物の基本的な物性記述子の計算。"""

from rdkit import Chem
from rdkit.Chem import Descriptors

from core.logging_utils import get_logger

logger = get_logger(__name__)


def calc_mol_weight(smiles: str) -> float | None:
    """SMILESから平均分子量(Da)を計算する。パースに失敗した場合はNoneを返す。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Failed to parse SMILES, skipping: %s", smiles)
        return None
    return Descriptors.MolWt(mol)
