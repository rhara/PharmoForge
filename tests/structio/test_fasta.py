from structio import parse_fasta


def test_parse_fasta_single_record(tmp_path):
    path = tmp_path / "single.fasta"
    path.write_text(">sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab\nMENFQK\nVEKI\n")

    records = parse_fasta(path)

    assert records == [("sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab", "MENFQKVEKI")]


def test_parse_fasta_multiple_records(tmp_path):
    path = tmp_path / "multi.fasta"
    path.write_text(">chain1\nGGGGAGGGGG\n>chain2\nMENFQKVEKI\n")

    records = parse_fasta(path)

    assert records == [("chain1", "GGGGAGGGGG"), ("chain2", "MENFQKVEKI")]
