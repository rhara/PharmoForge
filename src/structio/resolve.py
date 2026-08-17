"""複数の構造ファイルをコマンドライン引数から解決する(`--indir`と拡張子省略に対応)。

`click`の`UNPROCESSED`引数列(tokens)を受け取り、構造ファイルパスのリストに解決する。
複数の機能(`alignview`/`sequencealign`)で共通のCLI引数体系として使うため、
アトミックな技術要素として`structio`に切り出している。
"""

from pathlib import Path

import click

# 拡張子省略時に試す順(優先度順)
AUTO_EXTENSIONS = (".cif", ".mmcif", ".pdb")


def resolve_structure_tokens(
    tokens: tuple[str, ...], extensions: tuple[str, ...] = AUTO_EXTENSIONS
) -> list[Path]:
    """`--indir`と拡張子省略に対応しつつ、tokens列を構造ファイルパスのリストに解決する。

    `--indir DIR`は繰り返し指定可能で、以降のファイル名(拡張子省略可、既定では`.cif`優先、
    次に`.pdb`)を`DIR`配下から解決する。"/"を含む指定(または絶対パス)は`--indir`に
    よらずカレントディレクトリ相対 or 絶対パスとして扱う。

    `extensions`で拡張子省略時に試す拡張子と優先順位を呼び出し元ごとに変更できる
    (例: `sequencealign`は配列のみのFASTAも許容するため`.fasta`を追加する。`alignview`は
    PyMOLでの3次元表示が前提のため既定のまま構造ファイルの拡張子のみを使う)。
    """
    paths: list[Path] = []
    indir: Path | None = None
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "--indir":
            i += 1
            if i >= len(tokens):
                raise click.UsageError("--indir にはディレクトリを指定してください")
            indir = Path(tokens[i])
            i += 1
            continue
        paths.append(_resolve_one(token, indir, extensions))
        i += 1
    if not paths:
        raise click.UsageError("構造ファイルを1つ以上指定してください")
    return paths


def _resolve_one(token: str, indir: Path | None, extensions: tuple[str, ...]) -> Path:
    raw = Path(token)
    # "/"を含む(または絶対パスの)指定は--indirによらずそのまま(カレント相対 or 絶対)扱う
    if indir is None or "/" in token or raw.is_absolute():
        base = raw
    else:
        base = indir / raw

    if base.suffix:
        if not base.exists():
            raise click.UsageError(f"ファイルが見つかりません: {base}")
        return base

    for ext in extensions:
        candidate = base.with_suffix(ext)
        if candidate.exists():
            return candidate
    tried = ", ".join(str(base.with_suffix(ext)) for ext in extensions)
    raise click.UsageError(f"ファイルが見つかりません({tried})")
