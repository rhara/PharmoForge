"""UniProtエントリの取得と、創薬(構造生物学・メディシナルケミストリー)向け情報の抽出。"""

import requests

from core.logging_utils import get_logger

logger = get_logger(__name__)

UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def fetch_entry(accession: str) -> dict:
    """UniProt accession(例: P00533)の生エントリJSONを取得する。"""
    accession = accession.strip().upper()
    logger.info("Fetching UniProt entry for %s ...", accession)
    resp = requests.get(UNIPROT_ENTRY_URL.format(accession=accession), timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_fasta(accession: str) -> bytes:
    """UniProt accession(例: P00533)の正規配列をUniProt標準のFASTA形式で取得する。"""
    accession = accession.strip().upper()
    logger.info("Fetching UniProt FASTA for %s ...", accession)
    resp = requests.get(UNIPROT_FASTA_URL.format(accession=accession), timeout=30)
    resp.raise_for_status()
    return resp.content


def _location(feature: dict) -> tuple[int, int]:
    loc = feature["location"]
    return loc["start"]["value"], loc["end"]["value"]


def _features_by_type(entry: dict, feature_type: str) -> list[dict]:
    return [f for f in entry.get("features", []) if f.get("type") == feature_type]


def _first_comment_text(entry: dict, comment_type: str) -> str | None:
    for comment in entry.get("comments", []):
        if comment.get("commentType") == comment_type and comment.get("texts"):
            return comment["texts"][0]["value"]
    return None


def _xref_property(xref: dict, key: str) -> str | None:
    for prop in xref.get("properties", []):
        if prop.get("key") == key:
            return prop.get("value")
    return None


def extract_protein_info(entry: dict) -> dict:
    """UniProtの生エントリJSONから、創薬関連情報を平坦なdictに整理する。

    抽出項目: 配列・生物種・遺伝子/蛋白質名・機能概要・EC番号・活性部位/結合部位残基・
    翻訳後修飾(ジスルフィド結合・糖鎖化・修飾残基)・領域情報(膜貫通領域・シグナルペプチド・
    ドメイン)・キーワード・疾患関連情報・PDB構造(ID・測定手法・分解能)/AlphaFold DBの相互参照。
    """
    protein_desc = entry.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    protein_name = recommended.get("fullName", {}).get("value")
    ec_numbers = [e["value"] for e in recommended.get("ecNumbers", [])]

    genes = entry.get("genes", [])
    gene_name = genes[0].get("geneName", {}).get("value") if genes else None

    organism = entry.get("organism", {})
    sequence = entry.get("sequence", {})

    active_sites = [
        {"position": _location(f)[0], "description": f.get("description")}
        for f in _features_by_type(entry, "Active site")
    ]
    binding_sites = [
        {"start": _location(f)[0], "end": _location(f)[1], "ligand": f.get("ligand", {}).get("name")}
        for f in _features_by_type(entry, "Binding site")
    ]
    disulfide_bonds = [
        {"start": _location(f)[0], "end": _location(f)[1]}
        for f in _features_by_type(entry, "Disulfide bond")
    ]
    glycosylation_sites = [
        {"position": _location(f)[0], "description": f.get("description")}
        for f in _features_by_type(entry, "Glycosylation")
    ]
    modified_residues = [
        {"position": _location(f)[0], "description": f.get("description")}
        for f in _features_by_type(entry, "Modified residue")
    ]
    transmembrane_regions = [
        {"start": _location(f)[0], "end": _location(f)[1], "description": f.get("description")}
        for f in _features_by_type(entry, "Transmembrane")
    ]
    signal_peptides = [
        {"start": _location(f)[0], "end": _location(f)[1]} for f in _features_by_type(entry, "Signal")
    ]
    domains = [
        {"start": _location(f)[0], "end": _location(f)[1], "description": f.get("description")}
        for f in _features_by_type(entry, "Domain")
    ]

    diseases = [
        c["disease"]["diseaseId"]
        for c in entry.get("comments", [])
        if c.get("commentType") == "DISEASE" and c.get("disease")
    ]
    keywords = [k["name"] for k in entry.get("keywords", [])]

    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdb_structures = [
        {
            "id": x["id"],
            "method": _xref_property(x, "Method"),
            "resolution": _xref_property(x, "Resolution"),
        }
        for x in xrefs
        if x["database"] == "PDB"
    ]
    alphafold_ids = [x["id"] for x in xrefs if x["database"] == "AlphaFoldDB"]

    return {
        "accession": entry.get("primaryAccession"),
        "entry_name": entry.get("uniProtkbId"),
        "protein_name": protein_name,
        "gene_name": gene_name,
        "organism": organism.get("scientificName"),
        "taxon_id": organism.get("taxonId"),
        "sequence": sequence.get("value"),
        "length": sequence.get("length"),
        "mol_weight": sequence.get("molWeight"),
        "ec_numbers": ec_numbers,
        "function": _first_comment_text(entry, "FUNCTION"),
        "keywords": keywords,
        "diseases": diseases,
        "active_sites": active_sites,
        "binding_sites": binding_sites,
        "disulfide_bonds": disulfide_bonds,
        "glycosylation_sites": glycosylation_sites,
        "modified_residues": modified_residues,
        "transmembrane_regions": transmembrane_regions,
        "signal_peptide": signal_peptides[0] if signal_peptides else None,
        "domains": domains,
        "pdb_structures": pdb_structures,
        "alphafold_id": alphafold_ids[0] if alphafold_ids else None,
    }


def fetch_protein_info(accession: str) -> dict:
    """UniProt accessionから創薬関連情報を取得しdictで返す(fetch_entry + extract_protein_infoの組み合わせ)。"""
    return extract_protein_info(fetch_entry(accession))
