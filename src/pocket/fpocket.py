"""fpocketによる蛋白ポケット検出、および検出されたポケットに面する残基の抽出。"""

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from structio import parse_structure

from core.logging_utils import get_logger

logger = get_logger(__name__)

# fpocketの`<stem>_info.txt`内キー(正規化: strip+lower)→Pocketのフィールド名
_INFO_KEYS = {
    "score": "score",
    "druggability score": "druggability_score",
    "number of alpha spheres": "n_alpha_spheres",
    "volume": "volume",
    "total sasa": "total_sasa",
    "polar sasa": "polar_sasa",
    "apolar sasa": "apolar_sasa",
    "hydrophobicity score": "hydrophobicity_score",
}


@dataclass
class PocketResidue:
    """ポケットに面する残基1件(auth番号・チェーンID)。"""

    chain_id: str
    resnum: int
    resname: str


@dataclass
class Pocket:
    """fpocketが検出した1ポケット分の情報。"""

    pocket_id: int
    score: float
    druggability_score: float
    n_alpha_spheres: int
    volume: float
    total_sasa: float
    polar_sasa: float
    apolar_sasa: float
    hydrophobicity_score: float
    residues: list[PocketResidue] = field(default_factory=list)


def run_fpocket(structure_path: Path, work_dir: Path) -> list[Pocket]:
    """fpocketで`structure_path`(PDB/CIF)中のポケットを検出し、スコア降順で返す。

    fpocketは入力ファイルと同じディレクトリに出力(`<stem>_out/`)を書き出すコマンドラインツールのため、
    `structure_path`を`work_dir`にコピーしてから実行する。fpocketの生出力(PyMOL/VMD可視化スクリプト等)は
    `work_dir/<stem>_out/`にそのまま残る。ポケット周辺残基の抽出用に、fpocket自体には`-w pdb`
    (常にPDB形式で出力)を指定し、[`structio`](../structio)で読めるようにする(fpocket生成のmmCIFは
    `pdbx_PDB_model_num`列を欠くためProDyでパースできない)。
    """
    structure_path = Path(structure_path)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    stem = structure_path.stem
    input_path = work_dir / structure_path.name
    shutil.copyfile(structure_path, input_path)

    logger.info("Running fpocket on %s", structure_path)
    result = subprocess.run(
        ["fpocket", "-f", input_path.name, "-w", "pdb"],
        cwd=work_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"fpocket failed (exit {result.returncode}): {result.stderr.strip()}")

    out_dir = work_dir / f"{stem}_out"
    pocket_infos = _parse_info_file(out_dir / f"{stem}_info.txt")

    pockets = [
        Pocket(
            pocket_id=pocket_id,
            score=info["score"],
            druggability_score=info["druggability_score"],
            n_alpha_spheres=int(info["n_alpha_spheres"]),
            volume=info["volume"],
            total_sasa=info["total_sasa"],
            polar_sasa=info["polar_sasa"],
            apolar_sasa=info["apolar_sasa"],
            hydrophobicity_score=info["hydrophobicity_score"],
            residues=_extract_residues(out_dir / "pockets" / f"pocket{pocket_id}_atm.pdb"),
        )
        for pocket_id, info in pocket_infos.items()
    ]
    pockets.sort(key=lambda pocket: pocket.score, reverse=True)
    logger.info("Done: found %d pocket(s) in %s", len(pockets), structure_path)
    return pockets


def _parse_info_file(path: Path) -> dict[int, dict[str, float]]:
    """fpocketの`<stem>_info.txt`を、ポケットID→(正規化済みキー→値)のdictに変換する。

    行ごとの区切り記号(`:`)前後の空白の有無が項目によって揺れているため(例:
    `Score : \\t0.549`と`Hydrophobicity score:\\t46.429`が混在)、最初の`:`で分割して
    両辺をstripする方式で統一的に扱う。
    """
    infos: dict[int, dict[str, float]] = {}
    current: dict[str, float] | None = None
    for line in path.read_text().splitlines():
        header = re.match(r"Pocket (\d+) :", line.strip())
        if header:
            current = {}
            infos[int(header.group(1))] = current
            continue
        if current is None or ":" not in line:
            continue
        raw_key, _, raw_value = line.partition(":")
        key = _INFO_KEYS.get(raw_key.strip().lower())
        if key is not None:
            current[key] = float(raw_value.strip())
    return infos


def _extract_residues(atm_path: Path) -> list[PocketResidue]:
    """fpocketの`pocket<N>_atm.pdb`(ポケットに面する原子)から残基一覧(重複除去、chain/resnum昇順)を抽出する。"""
    atoms = parse_structure(atm_path)
    unique = {
        (str(chain_id), int(resnum)): str(resname)
        for chain_id, resnum, resname in zip(atoms.getChids(), atoms.getResnums(), atoms.getResnames())
    }
    return [
        PocketResidue(chain_id=chain_id, resnum=resnum, resname=resname)
        for (chain_id, resnum), resname in sorted(unique.items())
    ]
