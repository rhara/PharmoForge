import numpy as np

from structfit import fit_by_residue_number


def _atom_line(serial: int, resseq: int, chain: str, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:>5} {'CA':<4} {'GLY':>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {'C':>2}"
    )


def _write_pdb(path, coords, resseqs, chain="A"):
    lines = [_atom_line(i + 1, resseq, chain, *c) for i, (c, resseq) in enumerate(zip(coords, resseqs))]
    lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _random_coords(n, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-10, 10, size=(n, 3))


def _random_rotation(seed):
    rng = np.random.default_rng(seed)
    a, b, c = rng.uniform(0, 2 * np.pi, size=3)
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


def test_fit_by_residue_number_recovers_known_transform(tmp_path):
    target_coords = _random_coords(20, seed=1)
    resseqs = list(range(100, 120))

    rotation = _random_rotation(seed=2)
    translation = np.array([5.0, -3.0, 2.0])
    mobile_coords = (target_coords - translation) @ rotation  # rotation.T @ (target - t) の転置形

    target_path = tmp_path / "target.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_pdb(target_path, target_coords, resseqs)
    _write_pdb(mobile_path, mobile_coords, resseqs)

    result = fit_by_residue_number(mobile_path, target_path)

    assert result.n_residues == 20
    assert result.mobile_chain == "A"
    assert result.target_chain == "A"
    assert result.rmsd < 1e-3

    mobile_homogeneous = np.hstack([mobile_coords, np.ones((len(mobile_coords), 1))])
    fitted = (result.matrix @ mobile_homogeneous.T).T[:, :3]
    assert np.allclose(fitted, target_coords, atol=1e-2)


def test_fit_by_residue_number_uses_common_residues_only(tmp_path):
    coords = _random_coords(10, seed=3)
    target_path = tmp_path / "target.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_pdb(target_path, coords, list(range(1, 11)))
    _write_pdb(mobile_path, coords, list(range(6, 16)))  # resi 6-10のみ共通

    result = fit_by_residue_number(mobile_path, target_path)

    assert result.n_residues == 5


def test_fit_by_residue_number_raises_when_no_common_residues(tmp_path):
    coords = _random_coords(5, seed=4)
    target_path = tmp_path / "target.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_pdb(target_path, coords, list(range(1, 6)))
    _write_pdb(mobile_path, coords, list(range(100, 105)))  # 共通の残基番号なし

    try:
        fit_by_residue_number(mobile_path, target_path)
        assert False, "ValueErrorが発生しなかった"
    except ValueError:
        pass


def test_fit_by_residue_number_picks_more_complete_chain(tmp_path):
    coords_full = _random_coords(20, seed=5)
    resseqs_full = list(range(1, 21))

    target_path = tmp_path / "target.pdb"
    _write_pdb(target_path, coords_full, resseqs_full, chain="A")

    mobile_path = tmp_path / "mobile.pdb"
    lines = []
    # chain A: 5残基のみ(不完全なコピー)、chain B: 20残基すべて(完全なコピー)
    for i, (c, resseq) in enumerate(zip(coords_full[:5], resseqs_full[:5])):
        lines.append(_atom_line(i + 1, resseq, "A", *c))
    for i, (c, resseq) in enumerate(zip(coords_full, resseqs_full)):
        lines.append(_atom_line(100 + i, resseq, "B", *c))
    lines += ["TER", "END"]
    mobile_path.write_text("\n".join(lines) + "\n")

    result = fit_by_residue_number(mobile_path, target_path)

    assert result.mobile_chain == "B"
    assert result.n_residues == 20
