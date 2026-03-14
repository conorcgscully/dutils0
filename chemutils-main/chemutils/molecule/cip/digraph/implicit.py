from openeye import oechem

from .node import DigraphNode


class ImplicitHydrogenNode(DigraphNode):
    """Node representing an implicit hydrogen which does not have its own OpenEye atom."""

    @property
    def atomic_number(self) -> int:
        return 1

    @property
    def atomic_mass(self) -> float:
        return oechem.OEGetAverageWeight(1)  # type: ignore

    @property
    def name(self) -> str:
        return "H"
