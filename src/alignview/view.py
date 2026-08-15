"""複数のPDB/CIF構造をPyMOLで開き、先頭の構造に他を重ね合わせて表示する。"""

import shutil
import subprocess
import tempfile
from pathlib import Path

from core.logging_utils import get_logger
from structfit import fit_by_residue_number

logger = get_logger(__name__)

# align/superは`mobile, target`、cealignは`target, mobile`の順で引数を取る。
_ALIGN_COMMAND_TEMPLATES = {
    "align": "align {mobile}, {target}",
    "super": "super {mobile}, {target}",
    "cealign": "cealign {target}, {mobile}",
}

# 構造を見分けやすいよう、オブジェクトごとに順に割り当てる配色。
_OBJECT_COLORS = ["green", "cyan", "magenta", "yellow", "salmon", "skyblue", "orange", "purple"]

# `align`(配列ベース)専用: mobile構造がtargetの一部ドメインしか含まない場合、
# 相同性の高い別ドメイン(例: キナーゼ/偽キナーゼの対)に誤対応することがある。
# 同一蛋白の構造間では残基番号が(UniProt基準で)揃っている前提で、targetの探索範囲を
# mobileの残基番号レンジ±marginに絞り込み、誤対応を防ぐ。絞り込みで対応原子が
# 得られない場合(番号体系が揃っていない等)は絞り込みなしのalignにフォールバックする。
_RANGE_RESTRICTED_ALIGN_TEMPLATE = """\
python
try:
    _resis = set()
    cmd.iterate("{mobile} and polymer.protein and name CA", "_resis.add(int(resi))", space={{"_resis": _resis}})
    if not _resis:
        raise ValueError("no polymer residues in {mobile}")
    _lo, _hi = min(_resis) - {margin}, max(_resis) + {margin}
    _target_sel = "{target} and resi %d-%d" % (_lo, _hi)
    _result = cmd.align("{mobile}", _target_sel)
    if _result[1] == 0:
        raise ValueError("no atoms matched within restricted range")
    print("align-view: {mobile} -> {target} RMSD=%.3f (%d atoms, restricted to resi %d-%d)" % (_result[0], _result[1], _lo, _hi))
except Exception as _exc:
    print("align-view: restricted align failed for {mobile} (%s); falling back to unrestricted align" % _exc)
    _result = cmd.align("{mobile}", "{target}")
    print("align-view: {mobile} -> {target} RMSD=%.3f (%d atoms, unrestricted)" % (_result[0], _result[1]))
python end\
"""

# `number`専用: PyMOLのシーケンス/構造アルゴリズムを一切使わず、structfit(ProDy)で
# 残基番号の一致のみからCA原子を直接対応付けて求めた剛体変換を、そのままPyMOLに適用する。
# 最も直接的で決定的な対応付けだが、同一蛋白かつ構造間で残基番号(UniProt基準等)が
# 揃っていることが前提。
_NUMBER_TRANSFORM_TEMPLATE = """\
python
cmd.transform_object("{mobile}", {matrix}, transpose=0)
print("align-view: {mobile} -> {target} RMSD=%.3f (%d common residues, chain %s -> %s, number-based)" % ({rmsd}, {n_residues}, "{mobile_chain}", "{target_chain}"))
python end\
"""


def _unique_object_names(paths: list[Path]) -> list[str]:
    """ファイル名(拡張子なし)をPyMOLオブジェクト名にする。重複時は連番を付ける。"""
    names = []
    seen: dict[str, int] = {}
    for path in paths:
        stem = path.stem
        count = seen.get(stem, 0)
        seen[stem] = count + 1
        names.append(stem if count == 0 else f"{stem}_{count + 1}")
    return names


def build_pymol_script(paths: list[Path], method: str, align_margin: int = 20) -> str:
    """構造読み込み・配色・アラインを行うPyMOLスクリプト(.pml)本文を組み立てる。

    先頭のpathを基準(target)とし、残りをそれぞれ基準にアラインする。
    align_marginは`method="align"`の場合のみ使う(他ドメインへの誤対応防止、詳細は上記コメント参照)。
    `method="number"`の場合はここで(PyMOL起動前に)`structfit.fit_by_residue_number()`を呼び出し、
    実際に構造ファイルを読んで剛体変換を計算する(他のmethodと異なりファイルI/Oを伴う)。
    """
    names = _unique_object_names(paths)
    lines = [f"load {path.resolve()}, {name}" for path, name in zip(paths, names)]
    lines.append("show cartoon")

    colors = (_OBJECT_COLORS * (len(names) // len(_OBJECT_COLORS) + 1))[: len(names)]
    lines += [f"color {color}, {name}" for name, color in zip(names, colors)]

    target = names[0]
    if method == "align":
        lines += [
            _RANGE_RESTRICTED_ALIGN_TEMPLATE.format(mobile=mobile, target=target, margin=align_margin)
            for mobile in names[1:]
        ]
    elif method == "number":
        for path, mobile in zip(paths[1:], names[1:]):
            try:
                result = fit_by_residue_number(path, paths[0])
            except ValueError as exc:
                logger.warning("number-based fit failed for %s (%s); left unaligned", path, exc)
                lines.append(f'print("align-view: number-based fit failed for {mobile} ({exc})")')
                continue
            lines.append(
                _NUMBER_TRANSFORM_TEMPLATE.format(
                    mobile=mobile,
                    target=target,
                    matrix=result.matrix.flatten().tolist(),
                    rmsd=result.rmsd,
                    n_residues=result.n_residues,
                    mobile_chain=result.mobile_chain,
                    target_chain=result.target_chain,
                )
            )
    else:
        template = _ALIGN_COMMAND_TEMPLATES[method]
        lines += [template.format(mobile=mobile, target=target) for mobile in names[1:]]

    lines.append("zoom all")
    return "\n".join(lines) + "\n"


def launch_alignment_view(
    structure_paths: list[Path], method: str = "align", pymol_env: str = "pymol", align_margin: int = 20
) -> None:
    """指定した構造をPyMOLで開き、先頭の構造に他をアラインして表示する。

    rdkitのバージョン要件の都合でPyMOLは専用のconda/mamba環境(既定`pymol`)に
    インストールされている前提で、`mamba run -n <env> pymol`経由で起動する。
    """
    paths = [Path(p) for p in structure_paths]
    logger.info("Building PyMOL script for %d structure(s) (method=%s) ...", len(paths), method)
    script = build_pymol_script(paths, method, align_margin=align_margin)

    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as f:
        f.write(script)
        script_path = Path(f.name)

    mamba = shutil.which("mamba") or shutil.which("conda")
    if mamba is None:
        raise RuntimeError("mamba/condaコマンドが見つかりません。PATHを確認してください。")

    command = [mamba, "run", "-n", pymol_env, "pymol", str(script_path)]
    logger.info("Launching PyMOL (env=%s) ...", pymol_env)
    try:
        subprocess.run(command, check=True)
    finally:
        script_path.unlink(missing_ok=True)
    logger.info("PyMOL session closed.")
