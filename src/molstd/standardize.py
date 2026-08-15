"""ChEMBL Structure Pipeline (https://github.com/chembl/ChEMBL_Structure_Pipeline)
に倣った化合物構造の標準化。"""

from chembl_structure_pipeline import standardizer
from rdkit import Chem

from core.logging_utils import get_logger

logger = get_logger(__name__)


def standardize_smiles(smiles: str) -> str | None:
    """SMILESを標準化し、親構造(塩等を除いた形)のcanonical SMILESを返す。

    パースに失敗した場合はNoneを返す。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Failed to parse SMILES, skipping: %s", smiles)
        return None
    std_mol = standardizer.standardize_mol(mol)
    parent_mol, _ = standardizer.get_parent_mol(std_mol)
    return Chem.MolToSmiles(parent_mol)
