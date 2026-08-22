"""AutoDock Vinaによるドッキング実行(専用conda/mamba環境で起動)・結果パース。"""

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from meeko import gridbox

from core.logging_utils import get_logger

logger = get_logger(__name__)

_RESULT_RE = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)", re.MULTILINE)


@dataclass
class VinaPose:
    """出力PDBQT中の1ポーズ分のスコア(モード番号順、affinityが良い順に並ぶ)。"""

    mode: int
    affinity: float  # kcal/mol
    rmsd_lb: float
    rmsd_ub: float


@dataclass
class VinaResult:
    poses: list[VinaPose]
    output_path: Path

    @property
    def best_affinity(self) -> float:
        return self.poses[0].affinity


def calc_search_box(coords, padding: float = 4.0) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """座標配列(例: ポケット残基のCA座標)を包含するVina探索ボックスの中心・サイズを計算する。

    座標配列にポケット周辺残基のCAをそのまま渡す場合、それらは既にポケットの縁まで広がっているため、
    paddingは大きくしすぎない(既定4A)。大きくしすぎると探索空間が不必要に広がりドッキングが遅くなる。
    """
    return gridbox.calc_box(coords, padding)


def run_vina(
    rigid_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    output_path: Path,
    flex_pdbqt: Path | None = None,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    cpu: int | None = None,
    seed: int = 0,
    vina_env: str = "vina",
) -> VinaResult:
    """AutoDock Vinaでドッキングを実行し、出力ポーズのスコアを返す。

    vinaとrdkitはBoost.Pythonのビルドが競合するため`pharmoforge`本体環境には同居できず、
    vinaは専用のconda/mamba環境(既定`vina`)に置く前提で`mamba run -n <env> vina ...`経由で
    起動する(`pymolrun.run_pymol_script`と同じ回避パターン)。
    """
    mamba = shutil.which("mamba") or shutil.which("conda")
    if mamba is None:
        raise RuntimeError("mamba/condaコマンドが見つかりません。PATHを確認してください。")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        mamba, "run", "-n", vina_env, "vina",
        "--receptor", str(rigid_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--out", str(output_path),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
    ]
    if flex_pdbqt is not None:
        command += ["--flex", str(flex_pdbqt)]
    if cpu is not None:
        command += ["--cpu", str(cpu)]

    logger.info("Running Vina (env=%s): %s -> %s", vina_env, ligand_pdbqt, output_path)
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"vina failed (exit {result.returncode}): {result.stderr.strip()}")

    poses = parse_vina_output(output_path)
    logger.info(
        "Done: %s -> best affinity %.2f kcal/mol (%d pose(s))",
        ligand_pdbqt, poses[0].affinity, len(poses),
    )
    return VinaResult(poses=poses, output_path=output_path)


def parse_vina_output(output_path: Path) -> list[VinaPose]:
    """出力PDBQTの`REMARK VINA RESULT:`行からポーズごとのスコアを抽出する(モード番号順)。

    `run_vina`が内部で使うほか、既に実行済みの出力(キャッシュ)を再読み込みする際にも使う
    (ドッキングは1件あたり数分〜規模になりうるため、ノートブック側は既存ファイルをスキップし
    このパーサでスコアだけ読み直せるようにしている)。
    """
    text = Path(output_path).read_text()
    poses = []
    for i, match in enumerate(_RESULT_RE.finditer(text), start=1):
        affinity, rmsd_lb, rmsd_ub = (float(g) for g in match.groups())
        poses.append(VinaPose(mode=i, affinity=affinity, rmsd_lb=rmsd_lb, rmsd_ub=rmsd_ub))
    if not poses:
        raise ValueError(f"Vina出力からポーズを抽出できない: {output_path}")
    return poses
