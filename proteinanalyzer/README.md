# proteinanalyzer

UniProtから、計算化学・メディシナルケミストリーに役立つ蛋白情報を取得しJSONで保存する機能。

## 使い方

```bash
pf protein-info <識別子> --output <出力JSONファイル>
```

```bash
pf protein-info EGFR_HUMAN --output data/egfr_info.json
pf protein-info P00533 --output data/egfr_info.json
```

`<識別子>` にはUniProt entry name(例: `EGFR_HUMAN`)またはaccession(例: `P00533`)のいずれを与えてもよい
(entry nameは共通パッケージ[`src/idmap`](../src/idmap)でaccessionに変換される)。

### 出力

UniProt REST APIのエントリJSONから、以下の項目を抽出したJSONを出力する(`src/uniprot`の`extract_protein_info()`)。

| キー | 内容 |
| --- | --- |
| `accession` / `entry_name` | UniProt accession / entry name |
| `protein_name` / `gene_name` | 蛋白質名 / 遺伝子名 |
| `organism` / `taxon_id` | 生物種名 / NCBI Taxonomy ID |
| `sequence` / `length` / `mol_weight` | アミノ酸配列 / 配列長 / 分子量 |
| `ec_numbers` | EC番号(酵素の場合) |
| `function` | 機能概要(UniProtの`FUNCTION`コメント、最初のテキストのみ) |
| `keywords` | UniProtキーワード一覧 |
| `diseases` | 関連疾患名の一覧(`DISEASE`コメント) |
| `active_sites` | 活性部位残基(位置・説明) |
| `binding_sites` | リガンド結合部位残基(範囲・リガンド名) |
| `disulfide_bonds` | ジスルフィド結合(範囲) |
| `glycosylation_sites` | 糖鎖化部位(位置・説明) |
| `modified_residues` | その他の翻訳後修飾残基(位置・説明。リン酸化等) |
| `transmembrane_regions` | 膜貫通領域(範囲・説明) |
| `signal_peptide` | シグナルペプチド領域(範囲、なければ`null`) |
| `domains` | 機能ドメイン(範囲・説明) |
| `pdb_ids` | 相互参照されているPDB構造ID一覧 |
| `alphafold_id` | AlphaFold DBエントリID(なければ`null`) |

天然変異(natural variant)・変異導入実験(mutagenesis)等は現時点で対象外(必要になった時点で拡張)。

## 実装方針

- UniProtエントリの取得・情報抽出そのものは共通パッケージ[`src/uniprot`](../src/uniprot)で行う
  (`fetch_entry()`で生JSON取得、`extract_protein_info()`でdictに整理。アトミックな技術要素として独立)。
- 識別子(entry name/accession)の解決は共通パッケージ[`src/idmap`](../src/idmap)を利用する。
- 取得した蛋白情報をJSONファイルへ書き出す処理のみが本機能固有(`src/proteinanalyzer/report.py`)。

## テスト

```bash
pytest tests/proteinanalyzer tests/uniprot tests/idmap
```
