from seqextract import ChainSequence

from structfit import find_best_chain_for_residues

REFERENCE_SEQUENCE = "MDQYCILGRIGEGAHGIVFKAKHVETGEIVALKKVALRRLE"


def test_find_best_chain_for_residues_picks_chain_with_most_coverage():
    # chain A: 参照配列とほぼ同一(全残基が対応付く)。chain B: 無関係な短い配列(ほぼ対応しない)。
    chain_a = ChainSequence(chain_id="A", sequence=REFERENCE_SEQUENCE, resnums=list(range(1, len(REFERENCE_SEQUENCE) + 1)))
    chain_b = ChainSequence(chain_id="B", sequence="WWWWWWWWWW", resnums=list(range(1, 11)))

    result = find_best_chain_for_residues([chain_b, chain_a], REFERENCE_SEQUENCE, reference_resnums=[10, 11, 33])

    assert result is not None
    assert result.chain_id == "A"
    assert result.resnum_pairs == [(10, 10), (11, 11), (33, 33)]


def test_find_best_chain_for_residues_returns_none_when_no_chain_covers_anything():
    chain_b = ChainSequence(chain_id="B", sequence="WWWWWWWWWW", resnums=list(range(1, 11)))

    result = find_best_chain_for_residues([chain_b], REFERENCE_SEQUENCE, reference_resnums=[500, 501])

    assert result is None


def test_find_best_chain_for_residues_empty_mobile_chains():
    result = find_best_chain_for_residues([], REFERENCE_SEQUENCE, reference_resnums=[10])

    assert result is None
