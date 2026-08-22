"""相同蛋白の複数のX線構造にまたがる、共結晶化リガンドの接触残基のコンセンサスを求める。"""

from dataclasses import dataclass
from pathlib import Path

from seqalign import align_to_reference
from seqextract import get_chain_sequences
from structio import parse_structure

from core.logging_utils import get_logger

logger = get_logger(__name__)

# 結晶化添加物・修飾残基・イオン等、真のリガンドとみなさないヘテロ原子のresname。
DEFAULT_EXCLUDED_RESNAMES = frozenset({
    "DMS", "EDO", "GOL", "ACT", "ACE", "FMT", "PEG", "PG4", "MPD", "BME", "IPA",
    "SO4", "PO4", "NA", "CL", "MG", "CA", "ZN", "K", "MN", "NI", "CO",
    "CSD", "KCX", "OCS", "MSE", "SEP", "TPO", "PTR", "MLY", "ARG", "SER", "CYS",
})


@dataclass
class ConsensusLigandContacts:
    """複数構造にまたがるリガンド接触残基のコンセンサス集計結果(残基番号は`reference_sequence`基準)。"""

    anchor_resnums: list[int]
    n_ligands: int
    contact_counts: dict[int, int]
    min_count: int


def find_consensus_ligand_contacts(
    structure_paths: list[Path],
    reference_sequence: str,
    contact_distance: float = 4.5,
    min_ligand_atoms: int = 8,
    min_fraction: float = 0.2,
    excluded_resnames: frozenset[str] = DEFAULT_EXCLUDED_RESNAMES,
) -> ConsensusLigandContacts:
    """相同蛋白の複数のX線構造(`structure_paths`、PDB/CIF)に含まれる共結晶化リガンドの接触残基を、
    配列アラインメントで`reference_sequence`の番号にマッピングし、一定割合以上の構造で再現された
    残基だけをコンセンサスの接触残基として返す。

    共結晶化リガンドにはフラグメントスクリーニングのオフターゲットヒットや別ポケットに結合した
    ものも混ざりうるため、単一構造だけの接触では採用せず、`min_fraction`(リガンドを含む構造数に
    対する割合、既定0.2)以上で再現された残基だけを採用する(ただし最低2構造は要求する)。
    `min_ligand_atoms`未満のヘテロ残基は結晶化添加物由来のノイズとみなし無視する。
    `excluded_resnames`(既定`DEFAULT_EXCLUDED_RESNAMES`)に含まれるresnameは常にリガンド候補から除く。

    構造ごとに複数の鎖・複数のリガンドが含まれうる。各リガンドについて、`contact_distance`(Å)以内の
    蛋白CA原子を接触残基とし、そのリガンドを含む鎖の配列を`reference_sequence`にアラインメントして
    番号を変換する。同一構造内で複数のリガンドが同じ鎖に接触していれば、その構造は複数回カウントされる
    (1構造1カウントではなく1リガンド1カウント)。
    """
    contact_counts: dict[int, int] = {}
    n_ligands = 0
    logger.info("scanning %d structure(s) for co-crystallized ligands ...", len(structure_paths))
    for i, structure_path in enumerate(structure_paths, start=1):
        label = structure_path.stem
        atoms = parse_structure(structure_path)
        hetero = atoms.select("hetero and not water")
        ligand_resnames = (
            [r for r in sorted(set(hetero.getResnames())) if r not in excluded_resnames]
            if hetero is not None
            else []
        )

        found_any = False
        for resname in ligand_resnames:
            ligand_atoms = hetero.select(f"resname {resname}")
            if ligand_atoms is None or ligand_atoms.numAtoms() < min_ligand_atoms:
                continue
            contacts = atoms.select(
                f"name CA and same residue as (protein within {contact_distance} of "
                f"(hetero and not water and resname {resname}))"
            )
            if contacts is None:
                continue
            contact_resnums = {int(r) for r in contacts.getResnums()}
            chain_id = str(contacts.getChids()[0])
            chain = next((c for c in get_chain_sequences(atoms) if c.chain_id == chain_id), None)
            if chain is None:
                continue

            alignment = align_to_reference(reference_sequence, chain.sequence, chain.resnums)
            mapped = sorted({rp for rp, rn in alignment.query_resnum_by_ref_pos.items() if rn in contact_resnums})
            for resnum in mapped:
                contact_counts[resnum] = contact_counts.get(resnum, 0) + 1
            n_ligands += 1
            found_any = True
            logger.info(
                "[%d/%d] %s: ligand %s (%d atoms), %d contact residue(s) -> reference resnum %s",
                i, len(structure_paths), label, resname, ligand_atoms.numAtoms(), len(contact_resnums), mapped,
            )
        if not found_any:
            logger.info(
                "[%d/%d] %s: no ligand (or only excluded additives)", i, len(structure_paths), label,
            )

    min_count = max(2, round(min_fraction * n_ligands)) if n_ligands else 0
    anchor_resnums = sorted(r for r, c in contact_counts.items() if c >= min_count)
    logger.info(
        "consensus ligand contact residues (contacted in >= %d/%d ligand-bound structures): %s",
        min_count, n_ligands, anchor_resnums,
    )
    return ConsensusLigandContacts(
        anchor_resnums=anchor_resnums,
        n_ligands=n_ligands,
        contact_counts=contact_counts,
        min_count=min_count,
    )
