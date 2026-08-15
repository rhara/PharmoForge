"""RDKitの`MolStandardize`(rdMolStandardize)のみを用いた化合物構造の標準化。

外部パッケージ`chembl_structure_pipeline`には依存しない自前実装。
[ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)と同種の処理
(官能基・電荷の正規化、塩/溶媒和物の除去)をRDKit標準の`rdMolStandardize`部品で組み立てている。
"""

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from core.logging_utils import get_logger

logger = get_logger(__name__)

_uncharger = rdMolStandardize.Uncharger()


def standardize_smiles(smiles: str) -> str | None:
    """SMILESを標準化し、親構造(塩等を除いた形)のcanonical SMILESを返す。

    パースに失敗した場合はNoneを返す。

    処理手順:
      1. `rdMolStandardize.Cleanup`: サニタイズ・官能基/電荷の正規化・金属の解離。
      2. `rdMolStandardize.FragmentParent`: 共有結合の最大フラグメントを選び、塩・溶媒和物を除去。
      3. `rdMolStandardize.Uncharger`: 可能な範囲で電荷を中和する。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Failed to parse SMILES, skipping: %s", smiles)
        return None
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    mol = _uncharger.uncharge(mol)
    return Chem.MolToSmiles(mol)
