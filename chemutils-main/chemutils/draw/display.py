from typing import Any

from openeye import oechem

from .color import Color
from .draw import Format, draw_molecule
from .highlight import AtomBondPredicate, AtomPredicate, BondPredicate, Highlight
from .label import LabelFunction, get_enhanced_stereo_label

JUPYTER_WIDTH = 450
JUPYTER_HEIGHT = 150


def display_molecule(
    oemol: oechem.OEMolBase,
    /,
    *,
    width: int = JUPYTER_WIDTH,
    height: int = JUPYTER_HEIGHT,
    format: Format = "png",
    label: str | None = None,
    background_color: Color | None = "white",
    grayscale: bool = False,
    atom_labels: LabelFunction | None = get_enhanced_stereo_label,
    condense_abbreviations: bool = False,
    highlight: AtomBondPredicate | Highlight | list[Highlight] | None = None,
    highlight_atoms: AtomPredicate | None = None,
    highlight_bonds: BondPredicate | None = None,
) -> Any:
    """Display a molecule in Jupyter, using `draw_molecule`."""
    contents = draw_molecule(
        oemol,
        width=width,
        height=height,
        format=format,
        label=label,
        background_color=background_color,
        grayscale=grayscale,
        atom_labels=atom_labels,
        condense_abbreviations=condense_abbreviations,
        highlight=highlight,
        highlight_atoms=highlight_atoms,
        highlight_bonds=highlight_bonds,
    )
    match format.lower():
        case "png":
            # Avoid importing IPython if this function is unused
            from IPython.display import Image

            return Image(contents)  # type: ignore
        case "svg":
            # Avoid importing IPython if this function is unused
            from IPython.display import SVG

            return SVG(contents)  # type: ignore
        case _:
            raise ValueError(f"Failed to display molecule: invalid format: `{format}`.")
