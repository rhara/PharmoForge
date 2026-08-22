"""BLASTジョブのファイルキャッシュ・再開可能な実行(notebookでの繰り返し実行を想定)。"""

import pickle
from pathlib import Path

from core.logging_utils import get_logger

from .ncbi import fetch_hits, submit_blast, wait_for_blast

logger = get_logger(__name__)


def run_cached_blast(
    sequence: str,
    cache_dir: Path,
    program: str = "blastp",
    database: str = "pdb",
    entrez_query: str | None = None,
    poll_interval: float = 10.0,
    timeout: float = 600.0,
) -> list[dict]:
    """`blastsearch.blast_search`にファイルキャッシュと再開可能性を加えたもの(notebook向け)。

    BLASTジョブは数十秒〜数分かかることがあり、待機中にnotebookカーネルが再起動する等で
    処理が中断されても、ジョブ自体を再投函せずに再開できるようにする:

    - `cache_dir/blast_hits.pkl`が存在すればそれを返す(完了済みジョブの結果キャッシュ)。
    - なければ`cache_dir/blast_rid.txt`を確認し、既存のRequest ID(RID)があればジョブを
      再投函せず待機を再開する。なければ新規に投函しRIDを保存する。
    - 待機がタイムアウト(`TimeoutError`)した場合はRIDキャッシュを残したまま例外を伝播させる
      (再度呼べば同じRIDで待機を再開できる)。ジョブ自体が失敗(`RuntimeError`)した場合は
      RIDキャッシュを破棄する(再投函が必要なため)。
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    hits_cache = cache_dir / "blast_hits.pkl"
    rid_cache = cache_dir / "blast_rid.txt"

    if hits_cache.exists():
        logger.info("Using cached BLAST hits: %s", hits_cache)
        with open(hits_cache, "rb") as f:
            return pickle.load(f)

    if rid_cache.exists():
        rid = rid_cache.read_text().strip()
        logger.info("Resuming existing BLAST job %s (submitted previously) ...", rid)
    else:
        rid = submit_blast(sequence, program=program, database=database, entrez_query=entrez_query)
        rid_cache.write_text(rid)

    try:
        wait_for_blast(rid, poll_interval=poll_interval, timeout=timeout)
    except RuntimeError:
        # ジョブ自体が失敗した場合は再投函が必要なのでRIDキャッシュを破棄する
        rid_cache.unlink(missing_ok=True)
        raise
    # TimeoutErrorはRIDキャッシュを残したまま伝播させる(再度呼べば同じRIDで待機を再開できる)

    hits = fetch_hits(rid)
    with open(hits_cache, "wb") as f:
        pickle.dump(hits, f)
    rid_cache.unlink(missing_ok=True)
    logger.info("Saved BLAST hits to %s", hits_cache)
    return hits
