"""PyMOLスクリプトを専用conda/mamba環境で起動する。"""

import shutil
import subprocess
import tempfile
from pathlib import Path

from core.logging_utils import get_logger

logger = get_logger(__name__)


def run_pymol_script(script: str, pymol_env: str = "pymol") -> None:
    """PyMOLスクリプト文字列を一時ファイルに書き出し、GUIモードのPyMOLで実行する。

    rdkitのバージョン要件の都合でPyMOLは専用のconda/mamba環境(既定`pymol`)に
    インストールされている前提で、`mamba run -n <env> pymol <script>`経由で起動する。
    処理はPyMOLウィンドウを閉じるまでブロックする。実行後、一時ファイルは削除する。
    """
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
