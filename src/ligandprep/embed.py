"""SMILESからの3D配座生成、およびドッキング用PDBQT形式への変換。"""

from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from core.logging_utils import get_logger

logger = get_logger(__name__)

_EMBED_SEED = 0xF00D  # 再現性のため固定


def prepare_ligand_pdbqt(smiles: str, name: str, output_path: Path) -> Path:
    """SMILES1件から3D配座を生成し、ドッキング用PDBQTとしてoutput_pathに書き出す。

    互変異性体・電荷状態の標準化は呼び出し側の責務とする(`molstd.standardize_smiles`等)。
    ここでは3D配座生成(ETKDGv3)・MMFF94最適化・電荷/原子タイプ割当(meeko, Gasteiger)・
    PDBQT変換のみを行う。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILESをパースできない: {smiles}")
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", name)

    params = AllChem.ETKDGv3()
    params.randomSeed = _EMBED_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError(f"3D配座を生成できない: {name} ({smiles})")
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    molsetups = preparator.prepare(mol)
    pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(molsetups[0])
    if not success:
        raise ValueError(f"PDBQT変換に失敗した: {name}: {error_msg}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(pdbqt_string)
    logger.info("Prepared ligand PDBQT: %s -> %s", name, output_path)
    return output_path
