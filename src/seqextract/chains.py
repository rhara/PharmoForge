"""構造(ProDy Atomic)から蛋白チェーンごとの配列(1文字表記)を抽出する。

CA原子(観測された残基のみ)に基づく配列であり、電子密度が見えず欠損した
残基は含まれない(UniProtの完全配列とは異なりうる)。
"""

from dataclasses import dataclass

from prody.atomic.atomic import Atomic


@dataclass
class ChainSequence:
    """1蛋白チェーン分の配列情報。"""

    chain_id: str
    sequence: str  # 1文字表記、観測されたCA原子の順(=残基番号の昇順)
    resnums: list[int]  # sequence[i]に対応する残基番号

    @property
    def length(self) -> int:
        return len(self.sequence)


def get_chain_sequences(atoms: Atomic) -> list[ChainSequence]:
    """`atoms`に含まれる蛋白チェーンごとに配列を抽出する(チェーンID昇順)。

    水分子・リガンド等、蛋白でないチェーン(CA原子を持たないもの)は結果に含めない。
    """
    chains = []
    hv = atoms.getHierView()
    for chain in sorted(hv, key=lambda c: c.getChid()):
        ca = chain.select("protein and name CA")
        if ca is None:
            continue
        chains.append(
            ChainSequence(
                chain_id=chain.getChid(),
                sequence=ca.getSequence(),
                resnums=[int(r) for r in ca.getResnums()],
            )
        )
    return chains
