from molscaffold import compute_scaffold


def test_compute_scaffold_strips_side_chain():
    # トルエン: ベンゼン環スキャフォールド + メチル基は除去される
    assert compute_scaffold("Cc1ccccc1") == "c1ccccc1"


def test_compute_scaffold_keeps_ring_system():
    # ビフェニル(側鎖なし)はそのままスキャフォールド
    assert compute_scaffold("c1ccc(-c2ccccc2)cc1") == "c1ccc(-c2ccccc2)cc1"


def test_compute_scaffold_invalid_smiles_returns_none():
    assert compute_scaffold("not a smiles") is None
