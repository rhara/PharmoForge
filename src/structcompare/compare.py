"""構造(ProDy Atomic)間のチェーン単位配列比較。

ProDyの`matchChains`をラップする。`matchChains`はまず残基番号・残基名による
直接対応付け(`pwalign=False`)を試み、失敗時にBiopython(`Bio.Align.PairwiseAligner`、
未対応の場合は`Bio.pairwise2`)による配列アラインメントにフォールバックする
(`pwalign=True`)。

現行のprodyバージョン(2.6.1)では、`pwalign=True`のフォールバック時に返される
`AtomMap`の残基番号が実際の対応関係を反映しないことを実データで確認している
(例: AlphaFold全長モデルと結晶構造domain構築物を比較すると、重ね合わせ後の
RMSDが数十Åに達する=対応が破綻している)。返される%identity/%overlap自体は
`pwalign`の有無によらず一致する(検証済み)ため参考情報としては使えるが、
残基単位の対応(変異検出等)には使えない。また`pwalign=True`は入力の組み合わせに
よって内部で`StopIteration`を送出して丸ごと失敗することも実データで確認している
(prody側`getAlignedMatch`の既知の不具合と見られる)。

このため:
- `match_chains()`(%identity/%overlapの一覧)は`pwalign=True`で広く候補を探すが、
  上記の内部エラーで失敗した組は警告ログを出しつつ「対応なし」として扱う。
- `find_substitutions()`(残基単位の置換検出)は`pwalign=False`(残基番号ベース)
  のみを用い、対応が取れない場合は`matched=False`を返す(既存の`structfit`
  `--method number`と同じ前提: 同一蛋白の構造間ではPDBの残基番号が揃っている)。
"""

import warnings
from dataclasses import dataclass

from prody import matchChains
from prody.atomic.atomic import Atomic

from core.logging_utils import get_logger

logger = get_logger(__name__)

# matchChains の seqid/overlap は 0 を受け付けない(0 < x <= 100 が要件)ため、
# 実質的にフィルタなしとして機能する最小値を使う。
_MIN_SEQID = 1.0
_MIN_OVERLAP = 1.0


@dataclass
class ChainMatch:
    """2構造間で対応が取れたチェーンの組。"""

    chain_id_a: str
    chain_id_b: str
    seqid: float  # %
    overlap: float  # %
    n_matched: int  # 対応付けられた残基数


@dataclass
class ResidueSubstitution:
    """対応する位置で残基名が異なる箇所(アミノ酸置換)。"""

    resnum: int  # a側・b側で共通の残基番号(残基番号ベース対応のため一致)
    resname_a: str
    resname_b: str


@dataclass
class SubstitutionReport:
    """基準チェーンに対する残基単位の置換検出結果。"""

    matched: bool  # 残基番号ベースで対応するチェーンが見つかったか
    chain_id_b: str | None  # atoms_b側で対応が取れたチェーンID(matched=Falseの場合None)
    seqid: float | None  # %(matched=Falseの場合None)
    overlap: float | None  # %(matched=Falseの場合None)
    substitutions: list[ResidueSubstitution]


def _match(atoms_a: Atomic, atoms_b: Atomic, pwalign: bool):
    with warnings.catch_warnings():
        # pwalign=Trueの場合、このprodyバージョンではBio.Align.PairwiseAlignerの
        # 属性名不一致により常にBio.pairwise2へフォールバックし、非推奨警告が出る
        # (既知のprody側の挙動)。
        warnings.filterwarnings("ignore", message=r".*[Bb]io\.pairwise2.*")
        try:
            raw_matches = matchChains(atoms_a, atoms_b, seqid=_MIN_SEQID, overlap=_MIN_OVERLAP, pwalign=pwalign)
        except StopIteration:
            # prody.matchChains(pwalign=True)は入力の組み合わせによって内部の
            # getAlignedMatchでStopIterationを送出することがある(prody側の既知の不具合)。
            logger.warning("matchChains: prody内部エラーのため、このチェーンの組はスキップします")
            raw_matches = []
    # 蛋白でないチェーン同士(残基0件)が自明に100%一致として返ることがあるため除外する
    return [m for m in (raw_matches or []) if m[0].numAtoms() > 0]


def match_chains(atoms_a: Atomic, atoms_b: Atomic) -> list[ChainMatch]:
    """`atoms_a`と`atoms_b`のチェーンを全対全で比較し、対応が取れた組を
    %identity降順で返す(対応が全く取れなかった組は含まれない)。

    配列アラインメントによるフォールバックも許可する(`pwalign=True`)ため、
    残基番号体系が異なる構造間(例: 全長モデル対ドメインのみの結晶構造)でも
    %identity/%overlapを把握できる。
    """
    raw_matches = _match(atoms_a, atoms_b, pwalign=True)
    matches = [
        ChainMatch(
            chain_id_a=str(a.getChids()[0]),
            chain_id_b=str(b.getChids()[0]),
            seqid=float(seqid),
            overlap=float(overlap),
            n_matched=a.numAtoms(),
        )
        for a, b, seqid, overlap in raw_matches
    ]
    logger.info("matchChains: %d chain pair(s) matched", len(matches))
    return matches


def find_substitutions(atoms_a: Atomic, atoms_b: Atomic) -> SubstitutionReport:
    """`atoms_a`と`atoms_b`の最良一致チェーン同士で、残基番号が共通する位置の
    アミノ酸置換(残基名が異なる位置)を列挙する。

    残基番号ベースの対応付け(`pwalign=False`)のみを用いる(理由は本モジュール
    docstring参照)。対応するチェーンが見つからない場合(残基番号体系が異なる等)は
    `matched=False`を返す。
    """
    raw_matches = _match(atoms_a, atoms_b, pwalign=False)
    if not raw_matches:
        logger.info("find_substitutions: 残基番号ベースで対応するチェーンが見つかりませんでした")
        return SubstitutionReport(matched=False, chain_id_b=None, seqid=None, overlap=None, substitutions=[])

    a, b, seqid, overlap = raw_matches[0]  # seqid降順で先頭が最良一致
    substitutions = [
        ResidueSubstitution(resnum=int(resnum), resname_a=str(resname_a), resname_b=str(resname_b))
        for resnum, resname_a, resname_b in zip(a.getResnums(), a.getResnames(), b.getResnames())
        if resname_a != resname_b
    ]
    logger.info(
        "find_substitutions: %d substitution(s) found (seqid=%.1f%%, overlap=%.1f%%)",
        len(substitutions), seqid, overlap,
    )
    return SubstitutionReport(
        matched=True,
        chain_id_b=str(b.getChids()[0]),
        seqid=float(seqid),
        overlap=float(overlap),
        substitutions=substitutions,
    )
