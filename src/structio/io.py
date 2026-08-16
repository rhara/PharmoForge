"""PDB/CIF構造ファイルの読み書き(拡張子で自動判別、ProDyのAtomGroupを介する)。"""

from pathlib import Path

from prody.atomic.atomic import Atomic
from prody.proteins.ciffile import parseMMCIF, writeMMCIF
from prody.proteins.pdbfile import parsePDB, writePDB

_CIF_SUFFIXES = {".cif", ".mmcif"}


def parse_structure(path: Path) -> Atomic:
    """PDB/CIF形式の構造ファイルを拡張子で自動判別して読み込む。

    CIF形式は`unite_chains=True`で読み込み、チェーンIDに`auth_asym_id`
    (PyMOLやRCSBのWebサイト等で見えるチェーンID)を使う。ProDyの既定
    (`label_asym_id`)は同じauthチェーンに属する水分子・リガンド等を別チェーンとして
    細分化するため、PDB形式(auth相当のIDのみを持つ)との一貫性のためにも統一する。
    """
    path = Path(path)
    if path.suffix.lower() in _CIF_SUFFIXES:
        return parseMMCIF(str(path), unite_chains=True)
    return parsePDB(str(path))


def write_structure(atoms: Atomic, path: Path) -> None:
    """ProDyのAtomGroup(または選択結果)を、拡張子で自動判別したPDB/CIF形式で書き出す。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() in _CIF_SUFFIXES:
        writeMMCIF(str(path), atoms)
    else:
        writePDB(str(path), atoms)
