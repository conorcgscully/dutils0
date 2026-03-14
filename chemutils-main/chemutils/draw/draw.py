from pathlib import Path
from typing import Literal, overload

import numpy as np
from openeye import oechem
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import rdGeometry

from .color import Color, parse_color
from .highlight import (
    AtomBondPredicate,
    AtomPredicate,
    BondPredicate,
    Highlight,
)
from .label import LabelFunction, get_enhanced_stereo_label
from .to_rdkit import get_rdkit_molecule_and_highlights

Format = Literal["png", "svg"]


@overload
def draw_molecule(
    oemol: oechem.OEMolBase,
    /,
    *,
    width: int,
    height: int,
    format: Literal["png"],
    label: str | None = None,
    background_color: Color | None = "white",
    grayscale: bool = False,
    atom_labels: LabelFunction | None = None,
    condense_abbreviations: bool = False,
    highlight: AtomBondPredicate | Highlight | list[Highlight] | None = None,
    highlight_atoms: AtomPredicate | None = None,
    highlight_bonds: BondPredicate | None = None,
    aromaticity: Literal["kekulize", "dashed", "circle"] = "circle",
) -> bytes: ...


@overload
def draw_molecule(
    oemol: oechem.OEMolBase,
    /,
    *,
    width: int,
    height: int,
    format: Literal["svg"],
    label: str | None = None,
    background_color: Color | None = "white",
    grayscale: bool = False,
    atom_labels: LabelFunction | None = None,
    condense_abbreviations: bool = False,
    highlight: AtomBondPredicate | Highlight | list[Highlight] | None = None,
    highlight_atoms: AtomPredicate | None = None,
    highlight_bonds: BondPredicate | None = None,
    aromaticity: Literal["kekulize", "dashed", "circle"] = "circle",
) -> str: ...


def draw_molecule(
    oemol: oechem.OEMolBase,
    /,
    *,
    width: int = -1,
    height: int = -1,
    format: Format,
    label: str | None = None,
    background_color: Color | None = "white",
    grayscale: bool = False,
    atom_labels: LabelFunction | None = get_enhanced_stereo_label,
    condense_abbreviations: bool = False,
    highlight: AtomBondPredicate | Highlight | list[Highlight] | None = None,
    highlight_atoms: AtomPredicate | None = None,
    highlight_bonds: BondPredicate | None = None,
    aromaticity: Literal["kekulize", "dashed", "circle"] = "circle",
) -> str | bytes:
    """
    Draw a molecule to either a PNG or SVG.

    Args:
        oemol: OpenEye molecule to draw.
        label: Optional label to draw underneath the image.
        width: Image width in pixels.
        height: Image height in pixels.
        format: Either `"png"` or `"svg"`.
        background_color: Background color of the image. Defaults to white. Can be either a hex code or a color name.
            If None, the background will be transparent.
        grayscale: Should the image be grayscale?
        atom_labels: Optional function to generate a label for each atom.
        condense_abbreviations: Should common molecular fragments be condensed?
        highlight: Highlights to apply. Can be a single `highlight`, a list of `highlight`s, an `OEAtomBondSet`, or
            a list of atoms and bonds.
        highlight_atoms: Rule to select which atoms to highlight. Can be an atom, a list of atoms, an `OEAtomIter`,
            a `OEUnaryAtomBoolFunc`, a `OEAtomBondSet` or a `OEMatchBase`. For multiple or fancier highlights, use the
            `highlight` argument.
        highlight_bonds: Rule to select which bonds to highlight. Can be a bond, a list of bonds, an `OEBondIter`,
            a `OEUnaryBondBoolFunc`, a `OEAtomBondSet` or a `OEMatchBase`. For multiple or fancier highlights, use the
            `highlight` argument.
        aromaticity: How to draw aromatic bonds. Can be one of:
            - `"kekulize"`: Kekulize the molecule and draw aromatic bonds as single or double bonds.
            - `"dashed"`: Draw aromatic bonds as dashed lines.
            - `"circle"`: Draw aromatic rings as circles.

    Returns:
        PNG bytes if format is `png`, SVG string if format is `svg`.
    """
    mol, atom_highlights, bond_highlights = get_rdkit_molecule_and_highlights(
        oemol,
        condense_abbreviations=condense_abbreviations,
        atom_labels=atom_labels,
        highlight=highlight,
        highlight_atoms=highlight_atoms,
        highlight_bonds=highlight_bonds,
    )

    # Draw invisible highlights for atoms and bonds that are not highlighted
    # This means that images with and without highlights will be justified in the same way
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in atom_highlights:
            atom_highlights[atom.GetIdx()] = [(1, 1, 1, 0)]
    for bond in mol.GetBonds():
        if bond.GetIdx() not in bond_highlights:
            bond_highlights[bond.GetIdx()] = [(1, 1, 1, 0)]

    match format.lower():
        case "png":
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        case "svg":
            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        case _:
            raise ValueError(f"Failed to draw molecule: invalid format: `{format}`.")

    options = drawer.drawOptions()

    if background_color:
        options.setBackgroundColour(parse_color(background_color))
    else:
        options.clearBackground = False

    if grayscale:
        options.useBWAtomPalette()
    else:
        options.useCDKAtomPalette()

    options.scalingFactor = 40.0

    options.bondLineWidth = 2.8  # Make bonds thicker
    options.baseFontSize = 0.55  # Make font slightly smaller due to bold font
    options.fontFile = str(Path(__file__).parent / "NotoSans-SemiBold.ttf")
    options.additionalAtomLabelPadding = 0.05  # Add slight padding around atom labels

    # Recolor some atom types to be less harsh
    options.updateAtomPalette(
        {
            5: parse_color("#f99776"),
            7: parse_color("#2576e9"),
            8: parse_color("#ee2648"),
            9: parse_color("#99c70f"),
            15: parse_color("#ed7b11"),
            16: parse_color("#d2c404"),
            17: parse_color("#27ae1f"),
            34: parse_color("#eda20e"),
            35: parse_color("#a44e3b"),
            53: parse_color("#ae45c0"),
        }
    )

    options.unspecifiedStereoIsUnknown = True

    # Override default behaviour of kekulizing rings
    mol = rdMolDraw2D.PrepareMolForDrawing(
        mol,
        kekulize=aromaticity == "kekulize",
        addChiralHs=False,
        wedgeBonds=True,
        forceCoords=False,
        wavyBonds=False,
    )

    if aromaticity == "circle":
        circles = []

        for atom_ring, bond_ring in zip(
            mol.GetRingInfo().AtomRings(), mol.GetRingInfo().BondRings(), strict=False
        ):
            points = []
            is_aromatic = True
            for atom_idx in atom_ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                is_aromatic = is_aromatic and atom.GetIsAromatic()
                x, y, z = mol.GetConformer().GetAtomPosition(atom_idx)
                points.append((x, y))
            if is_aromatic:
                circles.append(get_circle_for_polygon(points))
                for bond_idx in bond_ring:
                    bond = mol.GetBondWithIdx(bond_idx)
                    bond.SetBondType(Chem.BondType.SINGLE)

    options.prepareMolsBeforeDrawing = False

    drawer.DrawMoleculeWithHighlights(
        mol,
        label if label else "",
        dict(atom_highlights),
        dict(bond_highlights),
        {atom.GetIdx(): 0.3 for atom in mol.GetAtoms()},
        {},
    )

    if aromaticity == "circle":
        drawer.SetFillPolys(False)
        for x, y, r in circles:
            r = r - 0.25  # How much to indent the rings by, in angstroms
            drawer.DrawEllipse(rdGeometry.Point2D(x - r, y - r), rdGeometry.Point2D(x + r, y + r))

    drawer.FinishDrawing()
    return drawer.GetDrawingText()  # type: ignore


def get_circle_for_polygon(points: list[tuple[float, float]], /) -> tuple[float, float, float]:
    # Note: this assumes that the points are a regular polygon
    # A more correct but more complicated solution is an inscribed circle
    points_np = np.array(points)
    midpoints = 0.5 * (points_np + np.roll(points_np, 1, axis=0))

    center = np.mean(points_np, axis=0)

    r = np.linalg.norm(midpoints - center, axis=1).min()

    return center[0], center[1], r
