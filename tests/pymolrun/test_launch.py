from pathlib import Path
from unittest.mock import patch

import pytest

from pymolrun.launch import run_pymol_script


def test_run_pymol_script_invokes_mamba_run_pymol():
    with (
        patch("pymolrun.launch.subprocess.run") as mock_run,
        patch("pymolrun.launch.shutil.which", return_value="/usr/bin/mamba"),
    ):
        run_pymol_script("load foo.pdb\n", pymol_env="pymol")

    args, kwargs = mock_run.call_args
    command = args[0]
    assert command[:4] == ["/usr/bin/mamba", "run", "-n", "pymol"]
    assert command[4] == "pymol"
    assert kwargs.get("check") is True

    script_path = Path(command[5])
    assert not script_path.exists()  # 実行後に一時ファイルが削除されている


def test_run_pymol_script_raises_when_mamba_not_found():
    with patch("pymolrun.launch.shutil.which", return_value=None):
        with pytest.raises(RuntimeError):
            run_pymol_script("load foo.pdb\n")
