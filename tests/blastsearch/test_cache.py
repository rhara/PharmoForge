from unittest.mock import patch

from blastsearch.cache import run_cached_blast


def test_run_cached_blast_uses_cached_hits_without_network(tmp_path):
    import pickle

    cache_dir = tmp_path
    (cache_dir / "blast_hits.pkl").write_bytes(pickle.dumps([{"subject_id": "sp|P11802|CDK4_HUMAN"}]))

    with patch("blastsearch.cache.submit_blast") as mock_submit:
        hits = run_cached_blast("MKV...", cache_dir)

    mock_submit.assert_not_called()
    assert hits == [{"subject_id": "sp|P11802|CDK4_HUMAN"}]


def test_run_cached_blast_submits_and_caches_on_first_run(tmp_path):
    with (
        patch("blastsearch.cache.submit_blast", return_value="RID123") as mock_submit,
        patch("blastsearch.cache.wait_for_blast") as mock_wait,
        patch("blastsearch.cache.fetch_hits", return_value=[{"subject_id": "sp|P11802|CDK4_HUMAN"}]) as mock_fetch,
    ):
        hits = run_cached_blast("MKV...", tmp_path, program="blastp", database="swissprot")

    mock_submit.assert_called_once()
    mock_wait.assert_called_once_with("RID123", poll_interval=10.0, timeout=600.0)
    mock_fetch.assert_called_once_with("RID123")
    assert hits == [{"subject_id": "sp|P11802|CDK4_HUMAN"}]
    assert (tmp_path / "blast_hits.pkl").exists()
    assert not (tmp_path / "blast_rid.txt").exists()  # 完了後はRIDキャッシュを削除する


def test_run_cached_blast_resumes_from_existing_rid_without_resubmitting(tmp_path):
    (tmp_path / "blast_rid.txt").write_text("RID999")

    with (
        patch("blastsearch.cache.submit_blast") as mock_submit,
        patch("blastsearch.cache.wait_for_blast") as mock_wait,
        patch("blastsearch.cache.fetch_hits", return_value=[]) as mock_fetch,
    ):
        run_cached_blast("MKV...", tmp_path)

    mock_submit.assert_not_called()
    mock_wait.assert_called_once_with("RID999", poll_interval=10.0, timeout=600.0)
    mock_fetch.assert_called_once_with("RID999")


def test_run_cached_blast_discards_rid_cache_on_job_failure(tmp_path):
    (tmp_path / "blast_rid.txt").write_text("RID999")

    with (
        patch("blastsearch.cache.submit_blast"),
        patch("blastsearch.cache.wait_for_blast", side_effect=RuntimeError("job failed")),
    ):
        try:
            run_cached_blast("MKV...", tmp_path)
        except RuntimeError:
            pass

    assert not (tmp_path / "blast_rid.txt").exists()


def test_run_cached_blast_keeps_rid_cache_on_timeout(tmp_path):
    (tmp_path / "blast_rid.txt").write_text("RID999")

    with (
        patch("blastsearch.cache.submit_blast"),
        patch("blastsearch.cache.wait_for_blast", side_effect=TimeoutError("still waiting")),
    ):
        try:
            run_cached_blast("MKV...", tmp_path)
        except TimeoutError:
            pass

    assert (tmp_path / "blast_rid.txt").read_text() == "RID999"
