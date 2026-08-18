import sqlite3

import pytest

from chembl import local


def _build_test_db(path):
    conn = sqlite3.connect(str(path))
    conn.executescript(
        """
        CREATE TABLE target_dictionary (tid INTEGER, target_type TEXT, pref_name TEXT, chembl_id TEXT);
        CREATE TABLE target_components (tid INTEGER, component_id INTEGER);
        CREATE TABLE component_sequences (component_id INTEGER, accession TEXT);
        CREATE TABLE activities (
            activity_id INTEGER, assay_id INTEGER, doc_id INTEGER, molregno INTEGER,
            standard_type TEXT, standard_value REAL, standard_units TEXT, pchembl_value REAL
        );
        CREATE TABLE assays (assay_id INTEGER, tid INTEGER, chembl_id TEXT);
        CREATE TABLE molecule_dictionary (molregno INTEGER, chembl_id TEXT, pref_name TEXT);
        CREATE TABLE compound_structures (molregno INTEGER, canonical_smiles TEXT);
        CREATE TABLE docs (doc_id INTEGER, chembl_id TEXT);

        INSERT INTO target_dictionary VALUES (1, 'SINGLE PROTEIN', 'CDK20', 'CHEMBL3559690');
        INSERT INTO target_dictionary VALUES (2, 'PROTEIN FAMILY', 'CDK family', 'CHEMBL3559691');
        INSERT INTO target_components VALUES (1, 100);
        INSERT INTO target_components VALUES (2, 100);
        INSERT INTO component_sequences VALUES (100, 'Q8IZL9');

        INSERT INTO activities VALUES (1, 10, 20, 1000, 'Kd', 8020.0, 'nM', 5.1);
        INSERT INTO activities VALUES (2, 10, 20, 1001, 'IC50', NULL, 'nM', NULL);

        INSERT INTO assays VALUES (10, 1, 'CHEMBL_ASSAY_1');
        INSERT INTO molecule_dictionary VALUES (1000, 'CHEMBL_MOL_1', 'Some Name');
        INSERT INTO molecule_dictionary VALUES (1001, 'CHEMBL_MOL_2', NULL);
        INSERT INTO compound_structures VALUES (1000, 'CCO');
        INSERT INTO docs VALUES (20, 'CHEMBL_DOC_1');
        """
    )
    conn.commit()
    conn.close()


def test_resolve_target_chembl_id_prefers_single_protein(tmp_path):
    db_path = tmp_path / "test.db"
    _build_test_db(db_path)

    assert local.resolve_target_chembl_id("Q8IZL9", db_path) == "CHEMBL3559690"


def test_resolve_target_chembl_id_raises_when_not_found(tmp_path):
    db_path = tmp_path / "test.db"
    _build_test_db(db_path)

    with pytest.raises(ValueError):
        local.resolve_target_chembl_id("NOTFOUND", db_path)


def test_fetch_activities_filters_pchembl_and_joins(tmp_path):
    db_path = tmp_path / "test.db"
    _build_test_db(db_path)

    records = local.fetch_activities("CHEMBL3559690", db_path)

    assert records == [
        {
            "molecule_chembl_id": "CHEMBL_MOL_1",
            "molecule_pref_name": "Some Name",
            "canonical_smiles": "CCO",
            "standard_type": "Kd",
            "standard_value": 8020.0,
            "standard_units": "nM",
            "pchembl_value": 5.1,
            "assay_chembl_id": "CHEMBL_ASSAY_1",
            "document_chembl_id": "CHEMBL_DOC_1",
        }
    ]


def test_fetch_activities_no_match_returns_empty(tmp_path):
    db_path = tmp_path / "test.db"
    _build_test_db(db_path)

    assert local.fetch_activities("CHEMBL_NOT_A_TARGET", db_path) == []
