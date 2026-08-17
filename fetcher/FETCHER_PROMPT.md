# fetcher 実装記録

このドキュメントは `fetcher` 機能を再現するための仕様記録。

## 目的

外部データベース(ChEMBL、RCSB PDB等)から、化合物・活性データや蛋白構造データを取得する。
PharmoForgeにおける「データ収集」機能を担う。

## CLI仕様

```
pf fetch <データ種別>=<識別子> --outdir <出力ディレクトリ>
```

- `<データ種別>` は現時点で `activity` (ChEMBL活性データ)、`structure` (構造データ、
  `,`区切りで複数可)をサポート。
- `--outdir` / `-o` は必須。出力先ディレクトリを指定する。ファイル名は識別子から自動的に決まる
  (都度明示的にファイル名を指定する必要はない)。
- `--type` (`cif`/`pdb`/`fasta`) は構造データのフォーマット指定に使う。省略時は`cif`。
- `--af` (フラグ、`structure=`専用): 取得元をRCSB PDBからAlphaFold DBに切り替える。
  `activity=`と併用するとエラー。

### activity: ChEMBL活性データ取得

```
pf fetch activity=CDK4_HUMAN --outdir data
```

- `<識別子>` にはUniProt entry name / UniProt accession / ChEMBL target idのいずれを与えてもよい。
  `src/idmap` の `resolve_chembl_target_id()` が自動判別して相互変換する。
- ChEMBL REST API (`https://www.ebi.ac.uk/chembl/api/data/activity.json`) を `target_chembl_id` と
  `pchembl_value__isnull=false` で絞り込み、`page_meta.next` を辿って全件取得する。
  `chembl_webresource_client` パッケージは使わず、`requests` による自前実装とする
  (`src/chembl/activity.py` の `fetch_activities()`。fetcher固有ではなくアトミックな技術要素として独立)。
- 取得したactivity recordの `canonical_smiles` をRDKit標準の`rdMolStandardize`のみで標準化する
  (`src/molstd` の `standardize_smiles()`、`rdMolStandardize.Cleanup()` → `FragmentParent()` → `Uncharger().uncharge()`
  で塩・電荷を正規化。外部パッケージ`chembl_structure_pipeline`には非依存の自前実装)。
- 標準化後のSMILES(`smiles`)が同一の化合物ごとに `pchembl_value` を集約し、
  median/mean/sd(sample stdev、1件のみの場合は空)、集約件数(`_n`)を計算する
  (`src/fetcher/activities.py` の `standardize_and_aggregate()`)。
  SMILESがパースできない・pChEMBL値が欠損している記録は集約対象から除外する。
- 出力行は`_median`の降順にソートする。
- 出力はTSV。列は `src/fetcher/activities.py` の `AGGREGATED_FIELDS` を参照
  (`smiles`, `_median`, `_mean`, `_sd`, `_n`)。
  元の`molecule_chembl_id`・`target_chembl_id`・アッセイ関連の列(`standard_type`等)は出力しない。
- 出力ファイル名は `<出力ディレクトリ>/<識別子>_activity.tsv`(上記の例では`data/CDK4_HUMAN_activity.tsv`)。

### structure: 構造データ取得(単体・複数、RCSB PDB / AlphaFold DB / UniProt)

```
pf fetch structure=9CSK --type=cif --outdir data
pf fetch structure=6P8F,7SJ3,9CSK --type pdb --outdir data
pf fetch structure=P61626 --af --type=cif --outdir data
pf fetch structure=CDK1_HUMAN --af --type=cif --outdir data/cdk1
pf fetch structure=P61626,CDK4_HUMAN --af --type pdb --outdir data
pf fetch structure=9CSK --type=fasta --outdir data
pf fetch structure=R1AB_SARS2 --type=fasta --outdir data
```

- `<識別子>` はカンマ区切りで複数指定できる(区切りがなければ単体扱い)。
- フォーマット(`cif`/`pdb`/`fasta`)は `--type` で指定する。省略時は`cif`。

`--type`が`cif`/`pdb`の場合:

- `--af`指定なし(既定): RCSB PDBのダウンロードエンドポイント
  `https://files.rcsb.org/download/{PDB_ID}.{fmt}` から取得する。`<識別子>`はPDB ID。
  出力ファイル名は `<出力ディレクトリ>/<識別子(大文字化)>.<fmt>`。
- `--af`指定あり: `<識別子>`にUniProt entry name(例: `CDK1_HUMAN`)またはaccession(例: `O14519`)を
  指定する。ダウンロード用accessionは`src/idmap`の`resolve_uniprot_accession()`が自動判別して解決する。
  AlphaFold DBのpredictionエンドポイント(`https://alphafold.ebi.ac.uk/api/prediction/{accession}`)から
  ダウンロードURL(`cifUrl`/`pdbUrl`)を都度解決する(バージョン番号をURLに固定しない)。
  出力ファイル名は識別子にentry nameを与えた場合でも常にUniProt accession(上記で解決したもの)
  に`_AF`を付けたものを使う: `<出力ディレクトリ>/<accession>_AF.<fmt>`
  (例: `pf fetch structure=CDK1_HUMAN --af --type=cif --outdir data/cdk1` →
  `data/cdk1/O14519_AF.cif`)。

`--type=fasta`の場合(詳細な経緯は[設計変更の経緯](#---typefastaの追加構造ファイルに加えfasta配列も取得できるように)を参照):

- RCSB/AlphaFold DBへは一切問い合わせず、常にUniProt本体
  (`https://rest.uniprot.org/uniprotkb/{accession}.fasta`、`src/uniprot/entry.py`の`fetch_fasta()`)
  から正規配列を直接取得する。`--af`は不要(付けても付けなくても結果は同じ)。
- `<識別子>`がPDB IDの形式であれば、そのPDBエントリに紐づくUniProt accessionをRCSB Data API
  (`data.rcsb.org`、`src/idmap/identifiers.py`の`pdb_id_to_uniprot_accessions()`)で解決する。
  1つのPDBエントリに複数の蛋白質(複合体)が含まれる場合は、accessionごとに1ファイル取得する。
- `<識別子>`がUniProt accessionの形式(`idmap.looks_like_uniprot_accession()`)、または`_`を含む
  (entry nameの形式、例: `R1AB_SARS2`)場合は、UniProt識別子とみなし`resolve_uniprot_accession()`
  で解決する。`--af`を明示的に付けた場合も常にこちらを使う。
  (PDB IDは4文字で`_`を含まず、UniProt accessionとも形式が異なるため、この判別に曖昧さはない)
- `--af`を指定した場合、fastaの取得元自体は変わらないため警告ログを出す(識別子解決は行われる)。
- 出力ファイル名は常に `<出力ディレクトリ>/<UniProt accession>.fasta`(`_AF`は付けない。RCSB/AlphaFold DB
  のどちらの経由でもなくUniProtから直接取得したものであることを表す)。
  複数識別子をまたいで同じaccessionが重複解決された場合は1回だけ取得する。

## 実装ファイル

- `src/fetcher/activities.py` — ChEMBL活性データの標準化+集約・TSV書き出し(fetcher固有の集計ロジック)
- `src/fetcher/cli.py` — `pf fetch` サブコマンド(データ種別・`--af`フラグ・`--type=fasta`のディスパッチ)
- `src/chembl/activity.py` — ChEMBL活性データの生取得(アトミックな技術要素として分離)
- `src/idmap/identifiers.py` — 蛋白識別子の相互マッピング(アトミックな技術要素として分離。`resolve_chembl_target_id()`は`activity=`、`resolve_uniprot_accession()`は`structure=`の`--af`指定時、`pdb_id_to_uniprot_accessions()`は`structure=`の`--type=fasta`(`--af`なし)時が使用。ダウンロード用accessionの解決と、出力ファイル名(`<accession>_AF.<fmt>`)にも同じaccessionを使う)
- `src/molstd/standardize.py` — 化合物構造の標準化(アトミックな技術要素として分離。`rdMolStandardize`のみの自前実装)
- `src/rcsb/download.py` — RCSB PDB構造データ取得(アトミックな技術要素として分離)
- `src/afdb/download.py` — AlphaFold DB構造データ取得(アトミックな技術要素として分離。`--type=fasta`時のみUniProt本体からの取得に委譲する)
- `src/uniprot/entry.py` — UniProtエントリ・配列の取得(アトミックな技術要素として分離。`fetch_fasta()`が`structure=` `--type=fasta`から使われる)

## 設計変更の経緯

### `structure-af=`/`structures-af=` → `structure=`/`structures=` + `--af`

当初はRCSB PDBとAlphaFold DBを別々のデータ種別(`structure`/`structures` vs
`structure-af`/`structures-af`)として実装していたが、ユーザーからの要望により、
データ種別は`structure=`/`structures=`に統一し、取得元をRCSB PDB(既定)/AlphaFold DB
(`--af`フラグ指定時)で切り替える方式に変更した。オプション体系(`--type`、複数形での
出力ディレクトリ扱い等)は共通のまま、`--af`の有無で内部的に`rcsb.fetch_structure`と
`afdb.fetch_structure`のどちらを呼ぶかを切り替える(`--af`指定時のみ識別子を
`idmap.resolve_uniprot_accession()`で解決する)。`--af`は`activities=`と併用するとエラー。

### `activities=`/`structure=`+`structures=` + `--output` → `activity=`/`structure=`(単複統合) + `--outdir`

その後さらにユーザーからの要望により、以下のように変更した。

- `activities=`(複数形)→`activity=`に変更。1標的の活性データをまとめて取得する処理内容は
  変わらないが、データ種別名を単数形に統一した。
- `structure=`(単体)/`structures=`(複数)を`structure=`に統合。識別子を`,`区切りで複数
  指定すれば複数取得、区切りがなければ単体取得と、識別子の書き方だけで単複を判別する。
- `--output`/`-o`(出力ファイルパスを都度明示、複数形の場合のみディレクトリ扱い)を廃止し、
  `--outdir`/`-o`(常に出力先ディレクトリを指定)に統一した。ファイル名は識別子から自動的に
  決まる: `activity=`は`<識別子>_activity.tsv`、`structure=`は`<識別子>.<fmt>`
  (`--af`指定時はaccessionを与えても常にUniProt entry nameを使う。`idmap.resolve_uniprot_entry_name()`
  で解決)。
- `--type`の省略時デフォルトを、旧来の「出力ファイル拡張子からの推定(既定cif)」から
  単純な既定値`cif`に変更した(出力ファイル名を都度指定しなくなったため、拡張子からの
  推定という概念自体が不要になった)。

### `--af`指定時の出力ファイル名: UniProt entry name → UniProt accession + `_AF`

上記の変更直後は`--af`指定時の出力ファイル名にUniProt entry name(`idmap.resolve_uniprot_entry_name()`
で解決)を使っていたが、ユーザーからの指摘により、UniProt accessionに`_AF`を付けたもの
(`<accession>_AF.<fmt>`)に変更した。これに伴い`resolve_uniprot_entry_name()`/
`accession_to_entry_name()`は不要になったため`src/idmap`から削除し、`--af`時の識別子解決は
`resolve_uniprot_accession()`(ダウンロード用accessionの解決と出力ファイル名の両方に使う)のみに戻した。

### バグ修正: `entry_name_to_accession()`があいまい検索により誤ったaccessionを返すことがある問題

`idmap.identifiers.entry_name_to_accession()`は当初UniProtの検索エンドポイント
(`GET /uniprotkb/search?query=id:<entry_name>`)を使い、先頭の検索結果(`results[0]`)を
採用していたが、この検索はトークン化された全文検索でありentry nameの完全一致ではないため、
無関係のエントリが先頭に返ることがあった。実例: `id:CDK1_HUMAN`で検索すると、完全一致する
`CDK1_HUMAN`(P06493、Cyclin-dependent kinase 1)より先に`CDKA1_HUMAN`(O14519、CDK2-associated
protein 1、全く別の蛋白)が返り、`structure=CDK1_HUMAN --af`が誤ってO14519の構造をダウンロード
していた。これはユーザーからの指摘で発覚した重大なバグ。

修正: UniProtのエントリ取得エンドポイント(`GET /uniprotkb/{id}.json`、`src/uniprot/entry.py`の
`UNIPROT_ENTRY_URL`と同種)はentry name/accessionのどちらを与えても該当accessionへ一意に解決される
(存在しない識別子の場合はHTTP 400)ため、こちらを使うように変更した。あいまいさのある検索
エンドポイントは`idmap`からは使わないようにした。

### `--type=fasta`の追加(構造ファイルに加えFASTA配列も取得できるように)

`--type`に`cif`/`pdb`に加え`fasta`を追加した。以下の3段階を経て現在の設計(fastaは常にUniProt本体
から直接取得し、`--af`はデータ取得元には影響しない)に至った。

**第1段階**: RCSB PDBは`files.rcsb.org/download`が`.fasta`拡張子に対応していないため、専用のfasta
エンドポイント(`https://www.rcsb.org/fasta/entry/{PDB_ID}`)に切り替えて取得するようにした
(`--af`指定時はAlphaFold DBのprediction APIレスポンスの`sequence`からFASTAを組み立てていた)。

**第2段階**: `--af`指定時のfastaについて、ユーザーからの指摘で以下の問題が発覚した。

- AlphaFold DBのモデルは全長ではなく断片(ドメイン単位)のことがある。実例:
  `structure=R1AB_SARS2 --af --type=fasta`(SARS-CoV-2のreplicase polyprotein 1ab、
  UniProt accession P0DTD1、全長7096残基)では、AlphaFold DBのprediction APIが返す
  `sequence`が126残基の断片モデルのものになっており、全長配列を取得できなかった。
- AlphaFold DBにモデルが存在しない識別子(`entries`が空)でも、UniProtエントリ自体は存在し
  配列取得が可能なケースがあり、AlphaFold API依存では不要にエラーになっていた。

これを受け、`--af`指定時のfastaはAlphaFold DBのprediction APIへは問い合わせず、UniProt本体の
fastaエンドポイント(`https://rest.uniprot.org/uniprotkb/{accession}.fasta`、`src/uniprot/entry.py`の
`fetch_fasta()`)から正規配列を直接取得する方式に変更した
(`afdb.download.fetch_structure`内で`fmt == "fasta"`を最初に分岐)。

**第3段階**: さらにユーザーから、「fastaはRCSB/AlphaFold DBという構造データベースとは無関係な
概念であり、UniProtから取得するのが確実」という指摘を受け、`--af`の有無によらずfastaは常に
UniProt本体から直接取得する設計に統一した。

- `--af`指定なし(PDB ID指定)の場合も、それまでのRCSBのfastaエンドポイントではなく、PDB IDに
  紐づくUniProt accessionを`idmap.pdb_id_to_uniprot_accessions()`(RCSB Data API `data.rcsb.org`で
  ポリマーエンティティごとのUniProt cross-referenceを解決)で求めてからUniProt本体を叩くように
  変更した。これに伴い`rcsb.download.fetch_structure`の`fasta`分岐は削除し、`cif`/`pdb`専用に戻した。
  1つのPDBエントリに複数の蛋白質(複合体)が含まれる場合は、accessionごとに複数ファイルを出力する
  (識別子1つに対しファイル1つという従来の前提が崩れるため、`fetcher/cli.py`側で`fmt == "fasta"`を
  別経路にした)。
- 出力ファイル名から`_AF`サフィックスを廃止し、常に`<UniProt accession>.fasta`にした
  (RCSB由来でもAlphaFold DB由来でもなく、UniProt本体から直接取得したものであるため)。
- `--af`はfastaのデータ取得元自体には影響しなくなった(常にUniProt)が、識別子の解決方法
  (PDB IDとして解決するかUniProt識別子として解決するか)には引き続き影響する。ユーザーの提案により、
  `--af --type=fasta`の組み合わせでは「取得元には影響しない」旨を警告ログとして出すようにした
  (識別子解決自体は行うため処理は継続する)。

**第4段階**: 第3段階の直後、`structure=R1AB_SARS2 --type=fasta`(`--af`なし)がエラーになる
(識別子`R1AB_SARS2`をPDB IDとして`data.rcsb.org`に問い合わせてしまい404になる)という指摘を受けた。
また、追加した警告ログ(`fetcher/cli.py`)が日本語だった点も指摘された
(ログ出力は`core.logging_utils`のロガーを使う箇所は英語、`click.UsageError`やdocstringは日本語、
という既存の使い分けに反していた)。以下の2点を修正した。

- `--af`の有無で識別子の解釈を切り替えるのをやめ、識別子の形式から自動判別するようにした:
  UniProt accessionの形式(`idmap.looks_like_uniprot_accession()`)に一致する、または`_`を含む
  (entry nameの形式、例: `R1AB_SARS2`)場合はUniProt識別子として`resolve_uniprot_accession()`で
  解決し、それ以外(PDB IDは4文字で`_`を含まず、UniProt accessionとも形式が異なる)はPDB IDとして
  `pdb_id_to_uniprot_accessions()`で解決する(`fetcher/cli.py`)。これにより`--af`は完全に不要になった
  (付けても付けなくても同じ結果になる。明示的に付けた場合は常にUniProt識別子として扱う経路を通す)。
- `fetcher/cli.py`に追加した警告ログを英語に変更した(`--af has no effect on FASTA output ...`)。
  同時に追加した`idmap.identifiers.pdb_id_to_uniprot_accessions()`内の`ValueError`メッセージも、
  同ファイルの既存の`ValueError`(英語)との整合性のため英語に変更した。

## テスト

ネットワークアクセスは `unittest.mock` でモックし、実際の外部APIへは接続しない。

```bash
pytest tests/fetcher tests/idmap tests/chembl tests/molstd tests/rcsb tests/afdb
```

## 動作例(サンプルデータ)

CDK4(ヒト、UniProt: P11802、ChEMBL: CHEMBL331)を題材にした活性データ取得例:

```bash
pf fetch activity=CDK4_HUMAN --outdir data
```
