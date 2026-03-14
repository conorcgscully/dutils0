from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING

from openeye import oechem

from .node import DigraphNode

if TYPE_CHECKING:
    from .digraph import HierarchicalDigraph


def make_subscript(value: int, /) -> str:
    """Convert a integer value to a subscripted string."""
    if value < 0:
        raise ValueError("Cannot convert negative number to subscript")
    return str(value).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))


class AtomNode(DigraphNode):
    """Node of a hierarchical digraph representing an atom."""

    atom: oechem.OEAtomBase
    """Underlying OpenEye atom."""

    def __init__(
        self,
        digraph: HierarchicalDigraph,
        atom: oechem.OEAtomBase,
    ):
        super().__init__(digraph=digraph)
        self.atom = atom

    @property
    def name(self) -> str:
        element = str(self.atom)[-2:].strip()
        if self.atom.GetAtomicNum() == 1:
            return element  # Don't show indices for H
        return f"{element}{make_subscript(self.atom.GetIdx())}"

    @property
    def parent(self) -> AtomNode | None:
        return super().parent  # type: ignore

    @property
    def ancestors(self) -> Collection[AtomNode]:
        return super().ancestors  # type: ignore

    @property
    def atomic_number(self) -> int:
        return self.atom.GetAtomicNum()  # type: ignore

    @property
    def atomic_mass(self) -> float:
        # OpenEye sets isotope to 0 when not present
        if (isotope := self.atom.GetIsotope()) > 0:
            return oechem.OEGetIsotopicWeight(self.atom.GetAtomicNum(), isotope)  # type: ignore
        # Atomic mass
        return oechem.OEGetAverageWeight(self.atomic_number)  # type: ignore

    @property
    def bond_to_parent(self) -> oechem.OEBondBase:
        if self.parent is None:
            raise ValueError("Atom Node has no parent")
        return self.atom.GetBond(self.parent.atom)

    @property
    def children(self) -> Collection[DigraphNode]:
        if len(super().children) == 0:
            self._explore_children()
        return super().children

    def _explore_children(self) -> None:
        """Generate the children of this node by exploring bonds and iterating implicit hydrogens."""
        for bond in self.atom.GetBonds():
            neighbour = bond.GetNbr(self.atom)
            # Create duplicates back to parent if multiple bond
            if self.parent is not None and neighbour == self.parent.atom:
                # IUPAC Blue Book, P-93.2.4
                # Treat all bonds from the root as single bonds (phosphates, etc.)
                if self.parent == self.digraph.root_node:
                    continue
                if bond.GetOrder() > 1:
                    for _ in range(bond.GetOrder() - 1):
                        self.digraph.add_duplicated(duplicated=self.parent, parent=self)
            # We've reached an atom previously encounted
            elif neighbour in [n.atom for n in self.ancestors]:
                self.digraph.add_duplicated(
                    duplicated=next(n for n in self.ancestors if n.atom == neighbour),
                    parent=self,
                )
                continue
            else:
                neighbour_node = self.digraph.add_atom(atom=neighbour, parent=self)
                # IUPAC Blue Book, P-93.2.4
                # Treat all bonds from the root as single bonds (phosphates, etc.)
                if self.digraph.root_node == self:
                    continue
                # Multiple bonds involve additional duplicated atoms
                if bond.GetOrder() > 1:
                    for _ in range(bond.GetOrder() - 1):
                        self.digraph.add_duplicated(duplicated=neighbour_node, parent=self)
        for _ in range(self.atom.GetImplicitHCount()):
            self.digraph.add_implicit_h(parent=self)
        for _ in range(3 - len(super().children)):
            self.digraph.add_phantom(parent=self)
