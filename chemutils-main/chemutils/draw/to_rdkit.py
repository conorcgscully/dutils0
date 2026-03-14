from openeye import oechem
from rdkit import Chem
from rdkit.Chem import rdAbbreviations, rdCoordGen, rdDepictor

from chemutils.rdmol_from_oemol import rdmol_from_oemol

from .highlight import (
    DEFAULT_HIGHLIGHT,
    AtomBondPredicate,
    AtomPredicate,
    BondPredicate,
    Highlight,
    get_highlights,
)
from .label import LabelFunction, get_enhanced_stereo_label


def get_rdkit_molecule_and_highlights(
    oemol: oechem.OEMolBase,
    /,
    *,
    condense_abbreviations: bool = False,
    atom_labels: LabelFunction | None = get_enhanced_stereo_label,
    highlight: AtomBondPredicate | Highlight | list[Highlight] | None = None,
    highlight_atoms: AtomPredicate | None = None,
    highlight_bonds: BondPredicate | None = None,
) -> Chem.Mol:
    oeatom_to_rdatom: dict[oechem.OEAtomBase, Chem.Atom] = {}
    oebond_to_rdbond: dict[oechem.OEBondBase, Chem.Bond] = {}
    mol = rdmol_from_oemol(oemol, atom_map=oeatom_to_rdatom, bond_map=oebond_to_rdbond)
    rdatom_to_oeatom = {v.GetIdx(): k for k, v in oeatom_to_rdatom.items()}

    if oemol.GetDimension() == 0:
        opts = rdCoordGen.CoordGenParams()
        opts.minimizerPrecision = opts.sketcherBestPrecision
        rdCoordGen.AddCoords(mol, opts)
        rdDepictor.StraightenDepiction(mol, minimizeRotation=True)

    if condense_abbreviations:
        mol = rdAbbreviations.CondenseMolAbbreviations(
            mol, rdAbbreviations.GetDefaultAbbreviations()
        )

    if atom_labels is not None:
        for atom in mol.GetAtoms():
            atom.SetProp("atomNote", str(atom_labels(rdatom_to_oeatom[atom.GetIdx()])))

    atom_highlights: dict[int, list[tuple[float, float, float]]] = {}
    bond_highlights: dict[int, list[tuple[float, float, float]]] = {}

    if highlight is not None and (highlight_atoms is not None or highlight_bonds is not None):
        raise ValueError(
            "Molecule highlighting error: can provide either `highlight_atoms` and/or`highlight_bonds`, or `highlight`."
        )

    if highlight is not None:
        if highlight_atoms is not None or highlight_bonds is not None:
            raise ValueError(
                "Molecule highlighting error: can provide either `highlight_atoms` and/or`highlight_bonds`, or `highlight`."
            )
        if not (
            isinstance(highlight, Highlight)
            or (isinstance(highlight, list) and all(isinstance(h, Highlight) for h in highlight))
        ):
            highlight = Highlight(atoms_and_bonds=highlight, color=DEFAULT_HIGHLIGHT)

    else:
        highlight = []
        if highlight_atoms:
            highlight.append(Highlight(atoms=highlight_atoms, color=DEFAULT_HIGHLIGHT))
        if highlight_bonds:
            highlight.append(Highlight(bonds=highlight_bonds, color=DEFAULT_HIGHLIGHT))

    if highlight:
        if isinstance(highlight, list):
            atom_highlights, bond_highlights = get_highlights(
                *highlight,
                rdmol=mol,
                oeatom_to_rdatom=oeatom_to_rdatom,
                oebond_to_rdbond=oebond_to_rdbond,
            )
        else:
            atom_highlights, bond_highlights = get_highlights(
                highlight,
                rdmol=mol,
                oeatom_to_rdatom=oeatom_to_rdatom,
                oebond_to_rdbond=oebond_to_rdbond,
            )

    return mol, atom_highlights, bond_highlights
