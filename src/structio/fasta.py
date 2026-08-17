"""FASTA形式の配列ファイルの読み込み。"""

from pathlib import Path

from Bio import SeqIO


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """FASTAファイルを(ヘッダー(先頭の">"を除く1行全体), 配列)のリストとして読み込む。"""
    return [(record.description, str(record.seq)) for record in SeqIO.parse(str(path), "fasta")]
