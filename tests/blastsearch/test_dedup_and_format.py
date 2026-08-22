from blastsearch import best_hit_per_accession, format_evalue


def test_format_evalue_zero():
    assert format_evalue(0.0) == "0"


def test_format_evalue_scientific_notation():
    assert format_evalue(9.47e-95) == "9.5e-95"


def test_best_hit_per_accession_keeps_lowest_evalue():
    hits = [
        {"subject_id": "sp|P11802|CDK4_HUMAN", "evalue": 1e-50},
        {"subject_id": "sp|P11802.1|CDK4_HUMAN", "evalue": 1e-60},  # 同accession、より良いevalue
        {"subject_id": "sp|P24941|CDK2_HUMAN", "evalue": 1e-40},
    ]

    result = best_hit_per_accession(hits)

    assert [h["accession"] for h in result] == ["P11802", "P24941"]
    assert result[0]["evalue"] == 1e-60


def test_best_hit_per_accession_excludes_given_accession():
    hits = [
        {"subject_id": "sp|Q8IZL9|CDK20_HUMAN", "evalue": 0.0},
        {"subject_id": "sp|P24941|CDK2_HUMAN", "evalue": 1e-40},
    ]

    result = best_hit_per_accession(hits, exclude_accession="Q8IZL9")

    assert [h["accession"] for h in result] == ["P24941"]
