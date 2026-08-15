from unittest.mock import MagicMock, patch

from uniprot import entry

SAMPLE_ENTRY = {
    "primaryAccession": "P00533",
    "uniProtkbId": "EGFR_HUMAN",
    "proteinDescription": {
        "recommendedName": {
            "fullName": {"value": "Epidermal growth factor receptor"},
            "ecNumbers": [{"value": "2.7.10.1"}],
        }
    },
    "genes": [{"geneName": {"value": "EGFR"}}],
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "sequence": {"value": "MRPSGTAGAA", "length": 10, "molWeight": 1111},
    "comments": [
        {"commentType": "FUNCTION", "texts": [{"value": "Receptor tyrosine kinase."}]},
        {
            "commentType": "DISEASE",
            "disease": {"diseaseId": "Lung cancer"},
        },
    ],
    "features": [
        {"type": "Active site", "location": {"start": {"value": 837}, "end": {"value": 837}}, "description": "Proton acceptor"},
        {
            "type": "Binding site",
            "location": {"start": {"value": 718}, "end": {"value": 726}},
            "ligand": {"name": "ATP"},
        },
        {"type": "Disulfide bond", "location": {"start": {"value": 31}, "end": {"value": 58}}},
        {"type": "Glycosylation", "location": {"start": {"value": 56}, "end": {"value": 56}}, "description": "N-linked"},
        {"type": "Modified residue", "location": {"start": {"value": 229}, "end": {"value": 229}}, "description": "Phosphoserine"},
        {
            "type": "Transmembrane",
            "location": {"start": {"value": 646}, "end": {"value": 668}},
            "description": "Helical",
        },
        {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 24}}},
        {
            "type": "Domain",
            "location": {"start": {"value": 712}, "end": {"value": 979}},
            "description": "Protein kinase",
        },
    ],
    "keywords": [{"name": "ATP-binding"}, {"name": "3D-structure"}],
    "uniProtKBCrossReferences": [
        {
            "database": "PDB",
            "id": "1IVO",
            "properties": [
                {"key": "Method", "value": "X-ray"},
                {"key": "Resolution", "value": "3.30 A"},
                {"key": "Chains", "value": "A/B=25-646"},
            ],
        },
        {
            "database": "PDB",
            "id": "2GS6",
            "properties": [{"key": "Method", "value": "NMR"}, {"key": "Chains", "value": "A=696-1022"}],
        },
        {"database": "AlphaFoldDB", "id": "P00533"},
        {"database": "STRING", "id": "9606.ENSP00000275493"},
    ],
}


def test_extract_protein_info_parses_all_fields():
    info = entry.extract_protein_info(SAMPLE_ENTRY)

    assert info["accession"] == "P00533"
    assert info["entry_name"] == "EGFR_HUMAN"
    assert info["protein_name"] == "Epidermal growth factor receptor"
    assert info["gene_name"] == "EGFR"
    assert info["organism"] == "Homo sapiens"
    assert info["taxon_id"] == 9606
    assert info["sequence"] == "MRPSGTAGAA"
    assert info["length"] == 10
    assert info["mol_weight"] == 1111
    assert info["ec_numbers"] == ["2.7.10.1"]
    assert info["function"] == "Receptor tyrosine kinase."
    assert info["keywords"] == ["ATP-binding", "3D-structure"]
    assert info["diseases"] == ["Lung cancer"]
    assert info["active_sites"] == [{"position": 837, "description": "Proton acceptor"}]
    assert info["binding_sites"] == [{"start": 718, "end": 726, "ligand": "ATP"}]
    assert info["disulfide_bonds"] == [{"start": 31, "end": 58}]
    assert info["glycosylation_sites"] == [{"position": 56, "description": "N-linked"}]
    assert info["modified_residues"] == [{"position": 229, "description": "Phosphoserine"}]
    assert info["transmembrane_regions"] == [{"start": 646, "end": 668, "description": "Helical"}]
    assert info["signal_peptide"] == {"start": 1, "end": 24}
    assert info["domains"] == [{"start": 712, "end": 979, "description": "Protein kinase"}]
    assert info["pdb_structures"] == [
        {"id": "1IVO", "method": "X-ray", "resolution": "3.30 A"},
        {"id": "2GS6", "method": "NMR", "resolution": None},
    ]
    assert info["alphafold_id"] == "P00533"


def test_extract_protein_info_handles_missing_optional_fields():
    minimal_entry = {
        "primaryAccession": "P61626",
        "uniProtkbId": "LYSC_HUMAN",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Lysozyme C"}}},
        "genes": [],
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "sequence": {"value": "MKALIV", "length": 6, "molWeight": 700},
        "comments": [],
        "features": [],
        "keywords": [],
        "uniProtKBCrossReferences": [],
    }

    info = entry.extract_protein_info(minimal_entry)

    assert info["gene_name"] is None
    assert info["ec_numbers"] == []
    assert info["function"] is None
    assert info["diseases"] == []
    assert info["signal_peptide"] is None
    assert info["alphafold_id"] is None
    assert info["pdb_structures"] == []


@patch("uniprot.entry.requests.get")
def test_fetch_entry_calls_expected_url(mock_get):
    mock_get.return_value = MagicMock(raise_for_status=lambda: None)
    mock_get.return_value.json.return_value = SAMPLE_ENTRY

    result = entry.fetch_entry("p00533")

    assert result == SAMPLE_ENTRY
    called_url = mock_get.call_args[0][0]
    assert called_url == "https://rest.uniprot.org/uniprotkb/P00533.json"


@patch("uniprot.entry.fetch_entry")
def test_fetch_protein_info_combines_fetch_and_extract(mock_fetch_entry):
    mock_fetch_entry.return_value = SAMPLE_ENTRY

    info = entry.fetch_protein_info("P00533")

    assert info["accession"] == "P00533"
    mock_fetch_entry.assert_called_once_with("P00533")
