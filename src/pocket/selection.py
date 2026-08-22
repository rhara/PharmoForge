"""fpocketが検出した複数のポケット候補から、既知のアンカー残基に基づいて1つを選ぶ。"""

from dataclasses import dataclass

from core.logging_utils import get_logger

from .fpocket import Pocket

logger = get_logger(__name__)


@dataclass
class PocketSelection:
    """アンカー残基との重なりに基づいて選ばれたポケット。"""

    pocket: Pocket
    overlap_resnums: list[int]

    @property
    def overlap(self) -> int:
        return len(self.overlap_resnums)


def select_pocket_by_anchor_overlap(
    pockets: list[Pocket], anchor_resnums: set[int], chain_id: str
) -> PocketSelection:
    """`anchor_resnums`(保存モチーフ・相同蛋白のリガンド接触残基等、生物学的根拠のある残基集合)との
    重なりが最大のポケットを選ぶ。

    fpocketのスコア(druggability等)は「目的のポケットかどうか」を直接表さない。特にリガンドを
    含まない予測構造では、他の表面ポケットがスコアで上回ることが普通に起こるため、スコアではなく
    既知のアンカー残基をどれだけ含むかでポケットを選ぶ。重なりが0件のポケットしかない場合は`ValueError`。
    """
    scored: list[tuple[int, Pocket, list[int]]] = []
    for p in pockets:
        resnums = {r.resnum for r in p.residues if r.chain_id == chain_id}
        overlap_resnums = sorted(resnums & anchor_resnums)
        scored.append((len(overlap_resnums), p, overlap_resnums))
        logger.info(
            "pocket %d: fpocket score=%.3f, anchor overlap=%d %s",
            p.pocket_id, p.score, len(overlap_resnums), overlap_resnums,
        )

    scored.sort(key=lambda t: t[0], reverse=True)
    best_overlap, best_pocket, best_overlap_resnums = scored[0]
    if best_overlap == 0:
        raise ValueError("アンカー残基を含むポケットが見つからない(fpocketの検出結果を確認してください)")

    logger.info(
        "Selected pocket %d (fpocket score=%.3f, anchor overlap=%d %s)",
        best_pocket.pocket_id, best_pocket.score, best_overlap, best_overlap_resnums,
    )
    return PocketSelection(pocket=best_pocket, overlap_resnums=best_overlap_resnums)
