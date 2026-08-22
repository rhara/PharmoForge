import pytest

from pocket import Pocket, PocketResidue
from pocket.selection import select_pocket_by_anchor_overlap


def _pocket(pocket_id, score, resnums, chain_id="A"):
    return Pocket(
        pocket_id=pocket_id,
        score=score,
        druggability_score=0.5,
        n_alpha_spheres=10,
        volume=100.0,
        total_sasa=50.0,
        polar_sasa=20.0,
        apolar_sasa=30.0,
        hydrophobicity_score=40.0,
        residues=[PocketResidue(chain_id=chain_id, resnum=r, resname="ALA") for r in resnums],
    )


def test_select_pocket_by_anchor_overlap_picks_best_overlap_not_best_score():
    # pocket 1はfpocketスコアが最も高いがアンカーと重ならない。pocket 2はスコアは低いが
    # アンカー2残基と重なる。fpocketスコアではなくアンカー重なりを優先すべき。
    pockets = [
        _pocket(1, score=0.9, resnums=[500, 501]),
        _pocket(2, score=0.1, resnums=[10, 11, 200]),
        _pocket(3, score=0.5, resnums=[10]),
    ]

    selection = select_pocket_by_anchor_overlap(pockets, anchor_resnums={10, 11, 33}, chain_id="A")

    assert selection.pocket.pocket_id == 2
    assert selection.overlap_resnums == [10, 11]
    assert selection.overlap == 2


def test_select_pocket_by_anchor_overlap_ignores_other_chains():
    pockets = [_pocket(1, score=0.5, resnums=[10, 11], chain_id="B")]

    with pytest.raises(ValueError):
        select_pocket_by_anchor_overlap(pockets, anchor_resnums={10, 11}, chain_id="A")


def test_select_pocket_by_anchor_overlap_raises_when_no_overlap():
    pockets = [_pocket(1, score=0.9, resnums=[500, 501])]

    with pytest.raises(ValueError):
        select_pocket_by_anchor_overlap(pockets, anchor_resnums={10, 11}, chain_id="A")
