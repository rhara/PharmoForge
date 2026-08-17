from unittest.mock import MagicMock, patch

import pytest

from blastsearch import ncbi


@patch("blastsearch.ncbi.requests.post")
def test_submit_blast_extracts_rid(mock_post):
    mock_post.return_value = MagicMock(
        text="    RID = ABC123XYZ\n    RTOE = 15\n", raise_for_status=lambda: None
    )

    rid = ncbi.submit_blast("MKV...", program="blastp", database="pdb")

    assert rid == "ABC123XYZ"
    called_kwargs = mock_post.call_args.kwargs
    assert called_kwargs["data"]["PROGRAM"] == "blastp"
    assert called_kwargs["data"]["DATABASE"] == "pdb"
    assert "ENTREZ_QUERY" not in called_kwargs["data"]


@patch("blastsearch.ncbi.requests.post")
def test_submit_blast_passes_entrez_query(mock_post):
    mock_post.return_value = MagicMock(text="    RID = ABC123XYZ\n", raise_for_status=lambda: None)

    ncbi.submit_blast("MKV...", database="swissprot", entrez_query="Homo sapiens[Organism]")

    called_kwargs = mock_post.call_args.kwargs
    assert called_kwargs["data"]["ENTREZ_QUERY"] == "Homo sapiens[Organism]"


@patch("blastsearch.ncbi.requests.post")
def test_submit_blast_raises_when_rid_missing(mock_post):
    mock_post.return_value = MagicMock(text="no rid here", raise_for_status=lambda: None)

    with pytest.raises(RuntimeError):
        ncbi.submit_blast("MKV...")


@patch("blastsearch.ncbi.time.sleep", return_value=None)
@patch("blastsearch.ncbi.requests.get")
def test_wait_for_blast_polls_until_ready(mock_get, _mock_sleep):
    mock_get.side_effect = [
        MagicMock(text="Status=WAITING\n", raise_for_status=lambda: None),
        MagicMock(text="Status=WAITING\n", raise_for_status=lambda: None),
        MagicMock(text="Status=READY\nThereAreHits=yes\n", raise_for_status=lambda: None),
    ]

    ncbi.wait_for_blast("ABC123XYZ", poll_interval=1.0, timeout=100.0)

    assert mock_get.call_count == 3


@patch("blastsearch.ncbi.requests.get")
def test_wait_for_blast_raises_on_failure_status(mock_get):
    mock_get.return_value = MagicMock(text="Status=FAILED\n", raise_for_status=lambda: None)

    with pytest.raises(RuntimeError):
        ncbi.wait_for_blast("ABC123XYZ", poll_interval=1.0, timeout=100.0)


@patch("blastsearch.ncbi.time.sleep", return_value=None)
@patch("blastsearch.ncbi.requests.get")
def test_wait_for_blast_raises_on_timeout(mock_get, _mock_sleep):
    mock_get.return_value = MagicMock(text="Status=WAITING\n", raise_for_status=lambda: None)

    with pytest.raises(TimeoutError):
        ncbi.wait_for_blast("ABC123XYZ", poll_interval=10.0, timeout=5.0)


def test_parse_pdb_subject_id_pipe_format():
    assert ncbi.parse_pdb_subject_id("pdb|6GZM|A") == ("6GZM", "A")


def test_parse_pdb_subject_id_underscore_format():
    assert ncbi.parse_pdb_subject_id("6gzm_A") == ("6GZM", "A")


def test_parse_uniprot_subject_id_plain_accession_with_version():
    assert ncbi.parse_uniprot_subject_id("Q8IZL9.1") == "Q8IZL9"


def test_parse_uniprot_subject_id_pipe_format():
    assert ncbi.parse_uniprot_subject_id("sp|P11802|CDK4_HUMAN") == "P11802"


@patch("blastsearch.ncbi.requests.get")
def test_fetch_hits_parses_tabular_output(mock_get):
    tabular = (
        "# BLASTP\n"
        "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, "
        "gap opens, q. start, q. end, s. start, s. end, evalue, bit score, % positives\n"
        "QUERY\tpdb|6GZM|A\t98.500\t250\t2\t0\t1\t250\t1\t250\t1e-150\t500.0\t99.0\n"
        "QUERY\tQ8IZL9.1\t85.000\t240\t10\t1\t1\t240\t5\t244\t1e-100\t400.0\t90.0\n"
    )
    mock_get.return_value = MagicMock(text=tabular, raise_for_status=lambda: None)

    hits = ncbi.fetch_hits("ABC123XYZ")

    assert hits == [
        {
            "subject_id": "pdb|6GZM|A",
            "identity": 98.5,
            "align_length": 250,
            "evalue": 1e-150,
            "bit_score": 500.0,
        },
        {
            "subject_id": "Q8IZL9.1",
            "identity": 85.0,
            "align_length": 240,
            "evalue": 1e-100,
            "bit_score": 400.0,
        },
    ]


@patch("blastsearch.ncbi.fetch_hits", return_value=[{"subject_id": "6GZM_A"}])
@patch("blastsearch.ncbi.wait_for_blast")
@patch("blastsearch.ncbi.submit_blast", return_value="ABC123XYZ")
def test_blast_search_combines_steps(mock_submit, mock_wait, mock_fetch):
    hits = ncbi.blast_search("MKV...", program="blastp", database="pdb")

    mock_submit.assert_called_once_with("MKV...", program="blastp", database="pdb", entrez_query=None)
    mock_wait.assert_called_once_with("ABC123XYZ", poll_interval=10.0, timeout=600.0)
    mock_fetch.assert_called_once_with("ABC123XYZ")
    assert hits == [{"subject_id": "6GZM_A"}]
