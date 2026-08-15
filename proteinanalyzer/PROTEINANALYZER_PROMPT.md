# proteinanalyzer 実装記録

このドキュメントは `proteinanalyzer` 機能を再現するための仕様記録。

## 目的

UniProtから、計算化学・メディシナルケミストリーに役立つ蛋白情報を取得しJSONで保存する。
PharmoForgeにおける「蛋白のシーケンス・構造解析」機能の第一歩。

## 抽出項目の選定方針

UniProtエントリの全項目のうち、創薬向け標準セットとして以下を選定した(ユーザーとの選択式確認で決定):

- 配列・生物種・遺伝子/蛋白質名・機能概要・EC番号
- 活性部位/結合部位残基
- 翻訳後修飾(ジスルフィド結合・糖鎖化・その他修飾残基)
- 領域情報(膜貫通領域・シグナルペプチド・ドメイン)
- キーワード・疾患関連情報
- PDB構造/AlphaFold DBの相互参照

天然変異(natural variant)・変異導入実験(mutagenesis)は情報量が多く用途も限定的なため、
現時点では対象外とした(必要になった時点で拡張)。

## CLI仕様

```
pf protein-info <識別子> --output <出力JSONファイル>
```

- `<識別子>` はUniProt entry name(例: `EGFR_HUMAN`)またはaccession(例: `P00533`)。
  `src/idmap`の`looks_like_uniprot_accession()`で判別し、entry nameの場合は`entry_name_to_accession()`で変換する。
- `--output` / `-o` は必須。

### 処理の流れ(`src/proteinanalyzer/cli.py`)

1. 識別子をaccessionに解決(`src/idmap`)。
2. `src/uniprot.fetch_protein_info(accession)` でUniProt REST API(`https://rest.uniprot.org/uniprotkb/{accession}.json`)
   からエントリJSONを取得し、創薬関連情報のdictに変換する。
3. `src/proteinanalyzer/report.write_protein_info_json()` でJSON(インデント付き、`ensure_ascii=False`)として書き出す。

## `src/uniprot` の実装詳細(`entry.py`)

- `fetch_entry(accession)`: UniProt REST APIから生エントリJSONを取得。
- `extract_protein_info(entry)`: 生JSONから以下を抽出する。
  - `proteinDescription.recommendedName` → `protein_name`、`ecNumbers` → `ec_numbers`
  - `genes[0].geneName` → `gene_name`
  - `organism.scientificName`/`taxonId`
  - `sequence.value`/`length`/`molWeight`
  - `comments`(`commentType == "FUNCTION"`の最初のテキスト → `function`、`commentType == "DISEASE"`の`disease.diseaseId`一覧 → `diseases`)
  - `features`(`type`別に`Active site`/`Binding site`/`Disulfide bond`/`Glycosylation`/`Modified residue`/`Transmembrane`/`Signal`/`Domain`を抽出。`location.start.value`/`end.value`を使用)
  - `keywords[].name` → `keywords`
  - `uniProtKBCrossReferences`(`database == "PDB"` → `pdb_ids`、`database == "AlphaFoldDB"`の最初のid → `alphafold_id`)
- `fetch_protein_info(accession)`: 上記2関数をまとめた入口。

## 実装ファイル

- `src/uniprot/entry.py` — UniProtエントリ取得・情報抽出(アトミックな技術要素として分離)
- `src/proteinanalyzer/report.py` — 蛋白情報dictのJSON書き出し(proteinanalyzer固有)
- `src/proteinanalyzer/cli.py` — `pf protein-info` サブコマンド

## テスト

```bash
pytest tests/proteinanalyzer tests/uniprot tests/idmap
```

`tests/uniprot/test_entry.py`は`extract_protein_info()`を実データ相当の合成JSON(EGFR由来の抜粋)で
実処理検証し、`fetch_entry()`のみ`unittest.mock`でネットワークをモックする。

## 動作例(サンプルデータ)

EGFR(ヒト、UniProt: P00533)を題材にした取得例(実際に動作確認済み。354件のPDB構造相互参照、
AlphaFold DBエントリ、活性部位・ATP結合部位・25件のジスルフィド結合等を正しく抽出できることを確認):

```bash
pf protein-info EGFR_HUMAN --output data/egfr_info.json
```
