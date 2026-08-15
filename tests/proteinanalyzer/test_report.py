import json

from proteinanalyzer.report import write_protein_info_json


def test_write_protein_info_json_writes_readable_json(tmp_path):
    info = {"accession": "P00533", "protein_name": "Epidermal growth factor receptor"}
    output = tmp_path / "out" / "egfr.json"

    write_protein_info_json(info, output)

    assert json.loads(output.read_text()) == info
