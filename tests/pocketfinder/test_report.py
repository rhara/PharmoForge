import io
import json

import pandas as pd

from pocket import Pocket, PocketResidue
from pocketfinder.report import (
    format_pockets_json,
    format_pockets_table,
    read_pockets_json,
    write_pockets_json,
    write_pockets_table,
)

_POCKET = Pocket(
    pocket_id=1,
    score=0.549,
    druggability_score=0.853,
    n_alpha_spheres=55,
    volume=441.369,
    total_sasa=64.374,
    polar_sasa=28.146,
    apolar_sasa=36.229,
    hydrophobicity_score=46.429,
    residues=[PocketResidue(chain_id="A", resnum=15, resname="TYR")],
)

_POCKET_2 = Pocket(
    pocket_id=2,
    score=0.255,
    druggability_score=0.001,
    n_alpha_spheres=17,
    volume=308.239,
    total_sasa=47.177,
    polar_sasa=29.062,
    apolar_sasa=18.114,
    hydrophobicity_score=41.429,
    residues=[
        PocketResidue(chain_id="A", resnum=3, resname="VAL"),
        PocketResidue(chain_id="B", resnum=9, resname="GLY"),
    ],
)


def test_format_pockets_json_structure():
    text = format_pockets_json("test.pdb", [_POCKET])
    data = json.loads(text)

    assert data["structure"] == "test.pdb"
    assert data["n_pockets"] == 1
    assert data["pockets"][0]["pocket_id"] == 1
    assert data["pockets"][0]["residues"] == [{"chain_id": "A", "resnum": 15, "resname": "TYR"}]


def test_write_pockets_json_creates_parent_dir(tmp_path):
    output = tmp_path / "nested" / "pockets.json"
    write_pockets_json("test.pdb", [_POCKET], output)

    data = json.loads(output.read_text())
    assert data["n_pockets"] == 1


def test_read_pockets_json_round_trips_write_pockets_json(tmp_path):
    output = tmp_path / "pockets.json"
    write_pockets_json("test.pdb", [_POCKET], output)

    pockets = read_pockets_json(output)

    assert pockets == [_POCKET]


def test_format_pockets_table_one_row_per_pocket_with_packed_residues():
    text = format_pockets_table([_POCKET, _POCKET_2])
    df = pd.read_csv(io.StringIO(text), sep="\t")

    assert list(df["pocket_id"]) == [1, 2]
    assert df.loc[0, "n_residues"] == 1
    assert df.loc[0, "residues"] == "A:15"
    assert df.loc[1, "n_residues"] == 2
    assert df.loc[1, "residues"] == "A:3,B:9"
    assert df.loc[0, "score"] == 0.549


def test_write_pockets_table_creates_parent_dir(tmp_path):
    output = tmp_path / "nested" / "pockets.tsv"
    write_pockets_table([_POCKET, _POCKET_2], output)

    df = pd.read_csv(output, sep="\t")
    assert len(df) == 2
    assert df.loc[1, "residues"] == "A:3,B:9"
