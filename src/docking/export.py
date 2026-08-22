"""ドッキング結果(Vina出力PDBQT)からの受容体コンフォメーション(PDB)・リガンドポーズ(SDF)の復元。

Vinaの`--out`出力自体にはリガンドポーズと可動側鎖(フレキシブル残基)の座標は含まれるが、受容体の
リジッド部分は含まれない(`--receptor`に渡した`_rigid.pdbqt`は不変のまま使い回される)。そのため、
インタラクション解析やMD初期構造として使える「その場でのポケット全体のコンフォメーション」を得るには、
リジッド部分の構造(`prepare_flexible_receptor`が書き出す`polymer_json`)と可動側鎖の座標(Vina出力)を
結合する必要がある。この結合をmeekoの`mk_export.py`(CLI)と同じ手順で行う。
"""

import copy
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from meeko import PDBQTMolecule, Polymer, RDKitMolCreate, export_pdb_updated_flexres

from core.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class ExportedPose:
    """1ポーズ分の書き出し結果(モード番号はVinaのスコア順、1始まり)。"""

    mode: int
    receptor_pdb: Path
    ligand_sdf: Path


def export_docked_poses(
    polymer_json: Path,
    vina_output_pdbqt: Path,
    output_dir: Path,
    name: str,
    modes: list[int] | None = None,
) -> list[ExportedPose]:
    """Vina出力(リガンド+フレキシブル受容体のPDBQT)から、モードごとに受容体のフルコンフォメーション
    (リジッド部分+可動側鎖、標準PDB形式)とリガンドポーズ(結合次数を復元したSDF)を書き出す。

    `polymer_json`は`prepare_flexible_receptor`が書き出した`<basename>.json`。`modes`省略時は
    全モードを書き出す(1始まり、Vinaのスコア順)。出力ファイルは
    `<output_dir>/<name>_mode<N>_receptor.pdb`・`<output_dir>/<name>_mode<N>_ligand.sdf`。
    """
    with open(polymer_json) as f:
        polymer = Polymer.from_json(f.read())

    pdbqt_mol = PDBQTMolecule.from_file(str(vina_output_pdbqt), skip_typing=True)
    sdf_string, failures = RDKitMolCreate.write_sd_string(pdbqt_mol, keep_flexres=False)
    if not sdf_string:
        raise ValueError(f"リガンドポーズをSDFに変換できない: {vina_output_pdbqt}")
    if failures:
        logger.warning("could not convert pose(s) %s to SDF for %s", sorted(failures), name)
    sdf_blocks = [block + "$$$$\n" for block in sdf_string.split("$$$$\n") if block.strip()]

    n_poses = pdbqt_mol._pose_data["n_poses"]
    if len(sdf_blocks) != n_poses:
        raise ValueError(
            f"SDFブロック数({len(sdf_blocks)})とポーズ数({n_poses})が一致しない: {vina_output_pdbqt}"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    exported = []
    for pose_id in range(n_poses):
        mode = pose_id + 1
        if modes is not None and mode not in modes:
            continue

        # export_pdb_updated_flexresはpolymer/pdbqt_molを書き換えるため、ポーズごとに複製する
        # (meekoのmk_export.py CLIと同じ手順)。
        single_pose = copy.deepcopy(pdbqt_mol)
        single_pose._positions = np.array([pdbqt_mol._positions[pose_id]])
        single_pose._pose_data["n_poses"] = 1
        single_pose._current_pose = 0
        receptor_pdb_string = export_pdb_updated_flexres(copy.deepcopy(polymer), single_pose)

        receptor_path = output_dir / f"{name}_mode{mode}_receptor.pdb"
        receptor_path.write_text(receptor_pdb_string)
        ligand_path = output_dir / f"{name}_mode{mode}_ligand.sdf"
        ligand_path.write_text(sdf_blocks[pose_id])

        exported.append(ExportedPose(mode=mode, receptor_pdb=receptor_path, ligand_sdf=ligand_path))

    logger.info("Exported %d pose(s) for %s -> %s", len(exported), name, output_dir)
    return exported
