"""ポケット一覧(list[Pocket])のJSON/TSV整形・ファイル書き出し。"""

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from pocket import Pocket, PocketResidue

from core.logging_utils import get_logger

logger = get_logger(__name__)


def _pockets_to_dict(structure_name: str, pockets: list[Pocket]) -> dict:
    return {
        "structure": structure_name,
        "n_pockets": len(pockets),
        "pockets": [asdict(pocket) for pocket in pockets],
    }


def format_pockets_json(structure_name: str, pockets: list[Pocket]) -> str:
    """ポケット一覧をJSON文字列(インデント付き、日本語等はエスケープしない)に整形する。"""
    return json.dumps(_pockets_to_dict(structure_name, pockets), ensure_ascii=False, indent=2)


def write_pockets_json(structure_name: str, pockets: list[Pocket], output: Path) -> None:
    """ポケット一覧をJSON(インデント付き、日本語等はエスケープしない)としてoutputに書き出す。"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(format_pockets_json(structure_name, pockets) + "\n")
    logger.info("Done: saved %d pocket(s) (%s) to %s", len(pockets), structure_name, output)


def read_pockets_json(path: Path) -> list[Pocket]:
    """`write_pockets_json`が書き出したJSON(`pf find-pocket`の出力)を`list[Pocket]`に復元する。"""
    data = json.loads(Path(path).read_text())
    return [
        Pocket(
            pocket_id=pocket["pocket_id"],
            score=pocket["score"],
            druggability_score=pocket["druggability_score"],
            n_alpha_spheres=pocket["n_alpha_spheres"],
            volume=pocket["volume"],
            total_sasa=pocket["total_sasa"],
            polar_sasa=pocket["polar_sasa"],
            apolar_sasa=pocket["apolar_sasa"],
            hydrophobicity_score=pocket["hydrophobicity_score"],
            residues=[PocketResidue(**residue) for residue in pocket["residues"]],
        )
        for pocket in data["pockets"]
    ]


def _residue_summary(residues: list[PocketResidue]) -> str:
    return ",".join(f"{residue.chain_id}:{residue.resnum}" for residue in residues)


def pockets_to_dataframe(pockets: list[Pocket]) -> pd.DataFrame:
    """ポケット一覧をDataFrame(1行1ポケット)に変換する。`residues`は`<chain>:<resnum>`をカンマ結合した1セルにまとめる。"""
    return pd.DataFrame(
        [
            {
                "pocket_id": pocket.pocket_id,
                "score": pocket.score,
                "druggability_score": pocket.druggability_score,
                "n_alpha_spheres": pocket.n_alpha_spheres,
                "volume": pocket.volume,
                "total_sasa": pocket.total_sasa,
                "polar_sasa": pocket.polar_sasa,
                "apolar_sasa": pocket.apolar_sasa,
                "hydrophobicity_score": pocket.hydrophobicity_score,
                "n_residues": len(pocket.residues),
                "residues": _residue_summary(pocket.residues),
            }
            for pocket in pockets
        ]
    )


def format_pockets_table(pockets: list[Pocket]) -> str:
    """ポケット一覧をTSV文字列(1行1ポケット)に整形する。"""
    return pockets_to_dataframe(pockets).to_csv(sep="\t", index=False)


def write_pockets_table(pockets: list[Pocket], output: Path) -> None:
    """ポケット一覧をTSV(1行1ポケット)としてoutputに書き出す。"""
    output.parent.mkdir(parents=True, exist_ok=True)
    pockets_to_dataframe(pockets).to_csv(output, sep="\t", index=False)
    logger.info("Done: saved %d pocket(s) as a table to %s", len(pockets), output)
