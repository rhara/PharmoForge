"""複数チェーンを持つ構造から、目的の残基集合を最もよくカバーするチェーンを選ぶ。"""

from dataclasses import dataclass

from seqalign import align_to_reference
from seqextract import ChainSequence

from core.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class BestChainCoverage:
    """`reference_resnums`を最もよくカバーするチェーンの選定結果。"""

    chain_id: str
    resnum_pairs: list[tuple[int, int]]  # (mobile側残基番号, reference側残基番号)


def find_best_chain_for_residues(
    mobile_chains: list[ChainSequence],
    reference_sequence: str,
    reference_resnums: list[int],
) -> BestChainCoverage | None:
    """`mobile_chains`(`seqextract.get_chain_sequences`の結果)の各チェーンを`reference_sequence`に
    アラインメントし、`reference_resnums`(基準配列側の残基番号の集合、例: ポケット周辺残基)を最も
    多くカバーするチェーンを選ぶ。

    結晶構造には非結晶学的対称(NCS)による複数コピーや、標的とは無関係な鎖(結晶化に使われた
    別蛋白等)が含まれることがあるため、目的の残基集合を最もよくカバーするチェーンを自動選択する
    (`structfit.fit_by_residue_pairs`にそのまま渡せる`resnum_pairs`を返す)。

    1件もカバーするチェーンがない場合は`None`。
    """
    best_chain_id: str | None = None
    best_pairs: list[tuple[int, int]] = []
    for chain in mobile_chains:
        alignment = align_to_reference(reference_sequence, chain.sequence, chain.resnums)
        pairs = [
            (alignment.query_resnum_by_ref_pos[r], r)
            for r in reference_resnums
            if r in alignment.query_resnum_by_ref_pos
        ]
        if len(pairs) > len(best_pairs):
            best_chain_id, best_pairs = chain.chain_id, pairs

    if best_chain_id is None:
        return None
    return BestChainCoverage(chain_id=best_chain_id, resnum_pairs=best_pairs)
