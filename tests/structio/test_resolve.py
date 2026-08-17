import click
import pytest

from structio.resolve import resolve_structure_tokens


def test_resolves_explicit_paths(tmp_path):
    p1 = tmp_path / "a.pdb"
    p2 = tmp_path / "b.cif"
    p1.write_text("ATOM\n")
    p2.write_text("data_b\n")

    assert resolve_structure_tokens((str(p1), str(p2))) == [p1, p2]


def test_indir_resolves_stems_with_auto_extension(tmp_path):
    (tmp_path / "a.cif").write_text("data_a\n")
    (tmp_path / "b.pdb").write_text("ATOM\n")

    result = resolve_structure_tokens(("--indir", str(tmp_path), "a", "b"))

    assert result == [tmp_path / "a.cif", tmp_path / "b.pdb"]


def test_indir_prefers_cif_over_pdb(tmp_path):
    (tmp_path / "a.cif").write_text("data_a\n")
    (tmp_path / "a.pdb").write_text("ATOM\n")

    result = resolve_structure_tokens(("--indir", str(tmp_path), "a"))

    assert result == [tmp_path / "a.cif"]


def test_indir_can_be_repeated(tmp_path):
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "a.cif").write_text("data_a\n")
    (dir2 / "b.pdb").write_text("ATOM\n")

    result = resolve_structure_tokens(("--indir", str(dir1), "a", "--indir", str(dir2), "b"))

    assert result == [dir1 / "a.cif", dir2 / "b.pdb"]


def test_slash_containing_token_ignores_indir(tmp_path):
    outside_dir = tmp_path / "outside"
    outside_dir.mkdir()
    outside_file = outside_dir / "c.cif"
    outside_file.write_text("data_c\n")
    indir = tmp_path / "indir"
    indir.mkdir()

    result = resolve_structure_tokens(("--indir", str(indir), str(outside_file)))

    assert result == [outside_file]


def test_missing_indir_value_raises():
    with pytest.raises(click.UsageError):
        resolve_structure_tokens(("--indir",))


def test_indir_stem_not_found_raises(tmp_path):
    with pytest.raises(click.UsageError):
        resolve_structure_tokens(("--indir", str(tmp_path), "missing"))


def test_no_tokens_raises():
    with pytest.raises(click.UsageError):
        resolve_structure_tokens(())


def test_default_extensions_do_not_resolve_fasta(tmp_path):
    (tmp_path / "a.fasta").write_text(">a\nMENFQKVEKI\n")

    with pytest.raises(click.UsageError):
        resolve_structure_tokens(("--indir", str(tmp_path), "a"))


def test_custom_extensions_can_include_fasta(tmp_path):
    (tmp_path / "a.fasta").write_text(">a\nMENFQKVEKI\n")

    result = resolve_structure_tokens(
        ("--indir", str(tmp_path), "a"), extensions=(".cif", ".pdb", ".fasta")
    )

    assert result == [tmp_path / "a.fasta"]


def test_custom_extensions_still_prefer_earlier_entries(tmp_path):
    (tmp_path / "a.cif").write_text("data_a\n")
    (tmp_path / "a.fasta").write_text(">a\nMENFQKVEKI\n")

    result = resolve_structure_tokens(
        ("--indir", str(tmp_path), "a"), extensions=(".cif", ".pdb", ".fasta")
    )

    assert result == [tmp_path / "a.cif"]


def test_explicit_fasta_extension_resolves_regardless_of_default_extensions(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    fasta_path.write_text(">a\nMENFQKVEKI\n")

    result = resolve_structure_tokens(("--indir", str(tmp_path), "a.fasta"))

    assert result == [fasta_path]
