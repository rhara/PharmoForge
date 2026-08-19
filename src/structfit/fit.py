"""同一蛋白の複数構造間で、残基番号の対応に基づく剛体重ね合わせ(rigid-body fit)を計算する。"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from prody import Transformation, calcRMSD, calcTransformation
from prody.atomic.atomic import Atomic

from core.logging_utils import get_logger
from structio import parse_structure

logger = get_logger(__name__)


@dataclass
class FitResult:
    """`target`に`mobile`を重ね合わせるための剛体変換の結果。"""

    matrix: np.ndarray  # 4x4行優先。v' = matrix @ [x, y, z, 1] (列ベクトル)
    rmsd: float
    n_residues: int
    mobile_chain: str
    target_chain: str


def _best_chain_pair(mobile, target) -> tuple[str, str, set[int]]:
    """共通の残基番号(CA原子)が最も多くなる(mobile鎖, target鎖)の組を選ぶ。

    結晶構造の非結晶学的対称(NCS)による複数コピー等、鎖ごとに構造の完全性が
    異なる場合でも、最も情報量の多い組み合わせを自動選択するため。
    """
    best: tuple[str, str, set[int]] | None = None
    for mobile_chain in sorted(set(mobile.getChids())):
        mobile_ca = mobile.select(f"name CA and protein and chain {mobile_chain}")
        if mobile_ca is None:
            continue
        mobile_resnums = {int(r) for r in mobile_ca.getResnums()}
        for target_chain in sorted(set(target.getChids())):
            target_ca = target.select(f"name CA and protein and chain {target_chain}")
            if target_ca is None:
                continue
            target_resnums = {int(r) for r in target_ca.getResnums()}
            common = mobile_resnums & target_resnums
            if best is None or len(common) > len(best[2]):
                best = (mobile_chain, target_chain, common)
    if best is None or not best[2]:
        raise ValueError("mobileとtargetの間に共通の残基番号(CA)が見つからない")
    return best


def fit_by_residue_number(mobile_path: Path, target_path: Path) -> FitResult:
    """残基番号が一致するCA原子同士を対応付け、mobileをtargetに重ね合わせる剛体変換を求める。

    同一蛋白の構造間ではPDBの残基番号が(UniProt基準等で)揃っている前提。
    配列アラインメントを一切行わないため対応が最も直接的で、番号が揃っている限り
    最も高精度な重ね合わせが得られる(反面、番号が揃っていない構造間には使えない)。
    """
    mobile = parse_structure(Path(mobile_path))
    target = parse_structure(Path(target_path))

    mobile_chain, target_chain, common_set = _best_chain_pair(mobile, target)
    common = sorted(common_set)

    mobile_ca = mobile.select(f"name CA and protein and chain {mobile_chain}")
    target_ca = target.select(f"name CA and protein and chain {target_chain}")
    mobile_index = {int(r): i for i, r in enumerate(mobile_ca.getResnums())}
    target_index = {int(r): i for i, r in enumerate(target_ca.getResnums())}

    mobile_coords = mobile_ca.getCoords()[[mobile_index[r] for r in common]]
    target_coords = target_ca.getCoords()[[target_index[r] for r in common]]

    transformation = calcTransformation(mobile_coords, target_coords)
    fitted = transformation.apply(mobile_coords.copy())
    rmsd = float(calcRMSD(fitted, target_coords))

    logger.info(
        "structfit: %s(chain %s) -> %s(chain %s): RMSD=%.3f (%d common residues)",
        mobile_path, mobile_chain, target_path, target_chain, rmsd, len(common),
    )
    return FitResult(
        matrix=transformation.getMatrix(),
        rmsd=rmsd,
        n_residues=len(common),
        mobile_chain=mobile_chain,
        target_chain=target_chain,
    )


def fit_by_residue_pairs(
    mobile_path: Path,
    target_path: Path,
    mobile_chain: str,
    target_chain: str,
    resnum_pairs: list[tuple[int, int]],
) -> FitResult:
    """`resnum_pairs`((mobile側残基番号, target側残基番号)の対応リスト)のCA原子同士を
    対応付け、mobileをtargetに重ね合わせる剛体変換を求める。

    `fit_by_residue_number`は同一蛋白の構造間(残基番号が一致する前提)専用だが、
    こちらは対応関係を呼び出し側が用意する(例: 別蛋白間の配列アラインメント結果から
    抽出した、ATP結合部位周辺residueのみの対応)ため、mobile/target間で残基番号体系が
    異なる構造間の重ね合わせにも使える。
    """
    mobile = parse_structure(Path(mobile_path))
    target = parse_structure(Path(target_path))

    mobile_ca = mobile.select(f"name CA and protein and chain {mobile_chain}")
    target_ca = target.select(f"name CA and protein and chain {target_chain}")
    if mobile_ca is None:
        raise ValueError(f"mobileにチェーン{mobile_chain}のCA原子が見つからない: {mobile_path}")
    if target_ca is None:
        raise ValueError(f"targetにチェーン{target_chain}のCA原子が見つからない: {target_path}")

    mobile_index = {int(r): i for i, r in enumerate(mobile_ca.getResnums())}
    target_index = {int(r): i for i, r in enumerate(target_ca.getResnums())}

    matched = [(m, t) for m, t in resnum_pairs if m in mobile_index and t in target_index]
    if not matched:
        raise ValueError(f"mobileとtargetの間に対応するCA原子の組が見つからない: {mobile_path} -> {target_path}")

    mobile_coords = mobile_ca.getCoords()[[mobile_index[m] for m, _ in matched]]
    target_coords = target_ca.getCoords()[[target_index[t] for _, t in matched]]

    transformation = calcTransformation(mobile_coords, target_coords)
    fitted = transformation.apply(mobile_coords.copy())
    rmsd = float(calcRMSD(fitted, target_coords))

    logger.info(
        "fit_by_residue_pairs: %s(chain %s) -> %s(chain %s): RMSD=%.3f (%d matched residue pairs)",
        mobile_path, mobile_chain, target_path, target_chain, rmsd, len(matched),
    )
    return FitResult(
        matrix=transformation.getMatrix(),
        rmsd=rmsd,
        n_residues=len(matched),
        mobile_chain=mobile_chain,
        target_chain=target_chain,
    )


def apply_fit(fit_result: FitResult, atoms: Atomic) -> Atomic:
    """`fit_result`(`fit_by_residue_number`/`fit_by_residue_pairs`の戻り値)の剛体変換を
    `atoms`(そのmobile構造由来のAtomic、選択結果でも可)に適用する。

    `Transformation.apply`はAtomPointer(選択結果)を渡した場合も、その親AtomGroup全体の
    座標を変換する(=チェーンID等のラベルは一切変更しない)ため、水を除いた選択結果を
    渡しても元の全チェーンの座標が正しく変換される。
    """
    return Transformation(fit_result.matrix).apply(atoms)
