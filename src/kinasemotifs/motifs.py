"""プロテインキナーゼドメインの保存モチーフを配列から検出する。"""

import re
from dataclasses import dataclass

from core.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class KinaseMotifs:
    """キナーゼドメインの保存モチーフの検出結果(1始まりの残基番号、区間は両端含む)。

    配列がこの正規表現に一致しない変種キナーゼ等では個々のモチーフがNoneになりうる
    (`find_kinase_motifs`は1つも見つからない場合にのみ`ValueError`にする)。
    """

    p_loop: tuple[int, int] | None  # GxGxxG
    catalytic_lys: int | None  # VAxK
    hrd: tuple[int, int] | None  # HRD
    dfg: tuple[int, int] | None  # DFG
    dfg_plus1: int | None  # DFGの直後(back pocketの壁を構成することが知られる位置)

    @property
    def anchor_resnums(self) -> set[int]:
        """検出できた全モチーフの残基番号をまとめた集合(ATP結合部位アンカーとして使う想定)。"""
        resnums: set[int] = set()
        if self.p_loop is not None:
            resnums |= set(range(self.p_loop[0], self.p_loop[1] + 1))
        if self.catalytic_lys is not None:
            resnums.add(self.catalytic_lys)
        if self.hrd is not None:
            resnums |= set(range(self.hrd[0], self.hrd[1] + 1))
        if self.dfg is not None:
            resnums |= set(range(self.dfg[0], self.dfg[1] + 1))
        if self.dfg_plus1 is not None:
            resnums.add(self.dfg_plus1)
        return resnums


def find_kinase_motifs(sequence: str) -> KinaseMotifs:
    """タンパク質配列からキナーゼドメインの保存モチーフを検出する。

    検出対象(いずれもATP結合部位を構成することが構造生物学的に確立している):
    P-loop(`G.G..G`、Gly-richループ、ATPのリン酸を上から覆う)、
    触媒Lys(`VA[LIV]K`、beta3鎖、ATPのalpha/betaリン酸と塩橋)、
    HRD(`HRD`、触媒ループ、基質のリン酸受容ヒドロキシル基を活性化)、
    DFG(`DFG`、AspがMg2+/ATPリン酸を配位)、
    DFG+1(DFGの直後の残基。モチーフ自体ではないが、多くのキナーゼでback pocketの壁を構成する
    ことが知られる位置のためアンカーに含める)。

    全て見つからない場合は`ValueError`。個々のモチーフが検出できないこと自体は許容する
    (`KinaseMotifs`の該当フィールドが`None`になる)。
    """
    p_loop = None
    m = re.search(r"G.G..G", sequence)
    if m:
        p_loop = (m.start() + 1, m.end())
        logger.info("P-loop (GxGxxG) motif: resnum %d-%d (%s)", p_loop[0], p_loop[1], m.group())

    catalytic_lys = None
    m = re.search(r"VA[LIV]K", sequence)
    if m:
        catalytic_lys = m.start() + 4
        logger.info("Catalytic Lys (VAxK) motif: resnum %d (%s)", catalytic_lys, m.group())

    hrd = None
    m = re.search(r"HRD", sequence)
    if m:
        hrd = (m.start() + 1, m.end())
        logger.info("Catalytic loop (HRD) motif: resnum %d-%d", hrd[0], hrd[1])

    dfg = None
    dfg_plus1 = None
    m = re.search(r"DFG", sequence)
    if m:
        dfg = (m.start() + 1, m.end())
        logger.info("DFG motif: resnum %d-%d", dfg[0], dfg[1])
        dfg_plus1 = m.end() + 1
        logger.info("DFG+1 (back pocket wall): resnum %d", dfg_plus1)

    motifs = KinaseMotifs(p_loop=p_loop, catalytic_lys=catalytic_lys, hrd=hrd, dfg=dfg, dfg_plus1=dfg_plus1)
    if not motifs.anchor_resnums:
        raise ValueError("キナーゼの保存モチーフ(P-loop/VAxK/HRD/DFG)が配列中に見つからない")
    return motifs
