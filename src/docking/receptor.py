"""ドッキング用受容体PDBQTの準備(指定残基の側鎖をフレキシブルにする)。"""

from dataclasses import dataclass
from pathlib import Path

import prody
from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer, ResidueChemTemplates

from core.logging_utils import get_logger

logger = get_logger(__name__)

_PRODY_PARSERS = {"pdb": prody.parsePDB, "cif": prody.parseMMCIF}


@dataclass
class FlexReceptor:
    """フレキシブル受容体のPDBQTファイル一式(Vinaの`--receptor`/`--flex`にそれぞれ対応)。"""

    rigid_pdbqt: Path
    flex_pdbqt: Path | None
    polymer_json: Path
    n_flexible_residues: int


def prepare_flexible_receptor(
    structure_path: Path,
    flexible_residues: list[tuple[str, int]],
    output_basename: Path,
) -> FlexReceptor:
    """構造ファイル(PDB/CIF)から受容体PDBQTを準備し、指定残基の側鎖をフレキシブルにする。

    `flexible_residues`((chain_id, resnum)のリスト)で指定した残基だけを可動側鎖として
    切り出し`<output_basename>_flex.pdbqt`に、それ以外は`<output_basename>_rigid.pdbqt`に
    書き出す。水素の付加・電荷割当(Gasteiger)はmeekoが内部で行う。

    受容体全体(リジッド+フレキシブル)のトポロジーは`<output_basename>.json`にも書き出す
    (`polymer.to_json()`)。Vinaの`--receptor`/`--flex`入力自体には座標更新のためリジッド部分の
    構造が含まれないため、ドッキング後に`export.export_docked_poses`で受容体のフルコンフォメーション
    (リジッド+ドッキング後の可動側鎖)をPDBとして復元する際に必要になる。
    """
    ext = Path(structure_path).suffix.lstrip(".").lower()
    parser = _PRODY_PARSERS.get(ext)
    if parser is None:
        raise ValueError(f"対応していない拡張子: {structure_path}")
    input_obj = parser(str(structure_path))

    templates = ResidueChemTemplates.create_from_defaults()
    mk_prep = MoleculePreparation()
    polymer = Polymer.from_prody(input_obj, templates, mk_prep, {}, [], False)

    for chain_id, resnum in flexible_residues:
        res_id = f"{chain_id}:{resnum}"
        if res_id not in polymer.monomers:
            raise ValueError(f"フレキシブル指定残基が受容体に見つからない: {res_id}")
        polymer.flexibilize_sidechain(res_id, mk_prep)

    rigid_pdbqt, flex_pdbqt_dict = PDBQTWriterLegacy.write_from_polymer(polymer)

    output_basename = Path(output_basename)
    output_basename.parent.mkdir(parents=True, exist_ok=True)
    rigid_path = output_basename.with_name(output_basename.name + "_rigid.pdbqt")
    rigid_path.write_text(rigid_pdbqt)

    flex_path = None
    if flex_pdbqt_dict:
        flex_path = output_basename.with_name(output_basename.name + "_flex.pdbqt")
        flex_path.write_text("".join(flex_pdbqt_dict.values()))

    json_path = output_basename.with_name(output_basename.name + ".json")
    json_path.write_text(polymer.to_json())

    logger.info(
        "Prepared flexible receptor: %s (%d flexible residue(s)) -> %s, %s, %s",
        structure_path, len(flexible_residues), rigid_path, flex_path, json_path,
    )
    return FlexReceptor(
        rigid_pdbqt=rigid_path,
        flex_pdbqt=flex_path,
        polymer_json=json_path,
        n_flexible_residues=len(flexible_residues),
    )
