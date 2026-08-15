"""スキャフォールド解析結果のレポート出力(TSV・グリッド画像・化合物レベルHTMLテーブル)。"""

import base64
import io
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from core.logging_utils import get_logger

logger = get_logger(__name__)

_COMPOUND_TABLE_STYLE = """
table { border-collapse: collapse; font-family: sans-serif; font-size: 13px; }
td, th { border: 1px solid #ccc; padding: 4px 8px; text-align: left; vertical-align: top; }
tr.scaffold-row td { background: #f0f0f0; font-weight: bold; }
tr.scaffold-row code { font-weight: normal; word-break: break-all; }
td.smiles { word-break: break-all; max-width: 320px; }
"""


def write_summary_tsv(summary: pd.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, sep="\t", index=False)
    logger.info("Wrote scaffold summary (%d scaffolds) to %s", len(summary), output)


def render_scaffold_grid(
    summary: pd.DataFrame,
    output: Path,
    top_n: int = 20,
    ascending: bool = False,
) -> None:
    """enrichmentの上位(ascending=False)または下位(ascending=True)top_n件のスキャフォールドを画像化する。"""
    rows = summary.sort_values("enrichment", ascending=ascending).head(top_n)

    mols = []
    legends = []
    for _, row in rows.iterrows():
        mol = Chem.MolFromSmiles(row["scaffold"])
        if mol is None:
            continue
        mols.append(mol)
        legends.append(
            f"n={int(row['n_total'])} high={int(row['n_high'])} low={int(row['n_low'])}\n"
            f"enrich={row['enrichment']:.3f} med={row['median_activity']:.2f}"
        )

    if not mols:
        logger.warning("No scaffolds to render, skipping %s", output)
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(250, 220), legends=legends, returnPNG=False)
    img.save(output)
    logger.info("Wrote scaffold grid image (%d scaffolds) to %s", len(mols), output)


def _mol_to_img_tag(smiles: str, size: tuple[int, int] = (180, 140)) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return f'<img src="data:image/png;base64,{b64}" width="{size[0]}" height="{size[1]}">'


def render_compound_table(
    df: pd.DataFrame,
    summary: pd.DataFrame,
    activity_col: str,
    output: Path,
    top_n: int = 20,
    ascending: bool = False,
    smiles_col: str = "smiles",
) -> None:
    """enrichment上位(ascending=False)/下位(ascending=True)top_n件のスキャフォールドについて、
    それに属する個々の化合物(置換基込み)を構造式付きのHTMLテーブルとして出力する。

    `df`は`scaffold`列・`bin`列・`smiles_col`列・`activity_col`列を持つこと
    (clustering.add_scaffolds/assign_activity_binsの出力)。
    スキャフォールドごとにグループ見出し行(構造式・n_total等)を挟み、
    グループ内は`activity_col`降順で化合物を並べる。
    """
    top_scaffolds = summary.sort_values("enrichment", ascending=ascending).head(top_n)

    row_htmls = []
    for _, srow in top_scaffolds.iterrows():
        scaffold = srow["scaffold"]
        scaffold_img = _mol_to_img_tag(scaffold, size=(160, 120))
        row_htmls.append(
            '<tr class="scaffold-row">'
            f"<td>{scaffold_img}</td>"
            f'<td colspan="3">'
            f"scaffold: <code>{scaffold}</code><br>"
            f"n_total={int(srow['n_total'])} n_high={int(srow['n_high'])} n_low={int(srow['n_low'])} "
            f"enrichment={srow['enrichment']:.3f}"
            "</td>"
            "</tr>"
        )
        compounds = df[df["scaffold"] == scaffold].sort_values(activity_col, ascending=False)
        for _, crow in compounds.iterrows():
            row_htmls.append(
                "<tr>"
                f"<td>{_mol_to_img_tag(crow[smiles_col])}</td>"
                f'<td class="smiles">{crow[smiles_col]}</td>'
                f"<td>{crow[activity_col]:.2f}</td>"
                f"<td>{crow['bin']}</td>"
                "</tr>"
            )

    output.parent.mkdir(parents=True, exist_ok=True)
    html = (
        "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>Scaffold compounds</title>"
        f"<style>{_COMPOUND_TABLE_STYLE}</style></head><body>\n"
        "<table>\n"
        f"<tr><th>structure</th><th>smiles</th><th>{activity_col}</th><th>bin</th></tr>\n"
        + "\n".join(row_htmls)
        + "\n</table>\n</body></html>\n"
    )
    output.write_text(html, encoding="utf-8")
    logger.info(
        "Wrote compound-level table (%d scaffolds, %d rows) to %s",
        len(top_scaffolds),
        len(row_htmls),
        output,
    )
