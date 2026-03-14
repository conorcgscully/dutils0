from __future__ import annotations

from collections.abc import Callable
from typing import Any, Generic, ParamSpec, TypeVar

from openeye import oechem
from rdkit import Chem

from chemutils.chemaxon.jvm import isjavainstance
from chemutils.chemaxon.jvm.chemaxon.struc import Molecule
from chemutils.chemaxon.molecule import cxmol_from_oemol
from chemutils.molecule import oemol_from_smiles
from chemutils.rdmol_from_oemol import rdmol_from_oemol

_T = TypeVar("_T", covariant=True)
_P = ParamSpec("_P")


class MolecularProperty(Generic[_T]):
    """
    Specific molecular property, which can be calculated using OpenEye or RDKit.

    Attributes:
        id: String ID used to identify this property as a column.
        oemol_func: Function to calculate this property given an OpenEye molecule.
        rdkit_func: Function to calculate this property given an RDKit molecule.

    """

    def __init__(
        self,
        id: str,
        oemol_func: Callable[[oechem.OEGraphMol], _T] | None = None,
        rdkit_func: Callable[[Chem.Mol], _T] | None = None,
        cxmol_func: Callable[[Molecule], _T] | None = None,
        composite_func: Callable[_P, _T] | None = None,
        composite_properties: dict[str, MolecularProperty[Any]] | None = None,
    ):
        self.id = id
        self.oemol_func = oemol_func
        self.rdkit_func = rdkit_func
        self.cxmol_func = cxmol_func
        self.composite_func = composite_func
        self.composite_properties = composite_properties
        # Check that only one is provided
        if (
            sum(
                [
                    oemol_func is not None,
                    rdkit_func is not None,
                    cxmol_func is not None,
                    composite_func is not None,
                ]
            )
            != 1
        ):
            raise ValueError(
                "Must only provide one of `oemol_func`, `rdkit_func`, `cxmol_func` or `composite_func`."
            )

    def __repr__(self) -> str:
        return f"<MolecularProperty '{self.id}'>"

    def get(self, mol: oechem.OEGraphMol | Chem.Mol | str, /) -> _T:
        """
        Get the value of this property for a molecule.

        The provided molecule can be an RDKit molecule, an OpenEye molecule or a SMILES string.

        Args:
            mol: Molecule, as either a `oechem.OEGraphMol`, `rdkit.Chem.Mol` or SMILES string.

        Returns:
            Calculated property value.

        Raises:
            ValueError: Cannot get property for object of unknown type.
        """
        if isinstance(mol, oechem.OEGraphMol):
            return self.get_openeye(mol)
        elif isinstance(mol, Chem.Mol):
            return self.get_rdkit(mol)
        elif isjavainstance(mol, Molecule):
            return self.get_chemaxon(mol)
        elif isinstance(mol, str):
            return self.get_openeye(oemol_from_smiles(mol))
        else:
            raise ValueError(f"Cannot get property for object of type {type(mol)}.")

    def get_composite(self, property_values: dict[MolecularProperty[Any], Any], /) -> _T:
        """
        Get the value of this property, given a dictionary of precomputed dependant properties.

        Args:
            property_values: Mapping of molecular properties this depends on to their values.

        Returns:
            Calculated property value.

        Raises:
            ValueError: Property is not a composite function.
        """
        if self.composite_func is None or self.composite_properties is None:
            raise ValueError("Property is not a composite function")
        return self.composite_func(  # type: ignore
            **{
                name: property_values[property]
                for name, property in self.composite_properties.items()
            }
        )

    def get_openeye(self, mol: oechem.OEGraphMol) -> _T:
        """
        Get the value of this property, given an OpenEye molecule.

        Args:
            mol: OpenEye molecule.

        Returns:
            Calculated property value.

        Raises:
            NotImplementedError: Property can not be calculated using OpenEye or RDKit.
        """
        if self.composite_func is not None and self.composite_properties is not None:
            return self.composite_func(  # type: ignore
                **{
                    name: property.get_openeye(mol)
                    for name, property in self.composite_properties.items()
                }
            )
        if self.cxmol_func is not None:
            return self.get_chemaxon(cxmol_from_oemol(mol))
        elif self.rdkit_func is not None:
            return self.get_rdkit(rdmol_from_oemol(mol))
        elif self.oemol_func is not None:
            return self.oemol_func(mol)
        else:
            raise NotImplementedError(
                "Property is not calculated using RDKit, Chemaxon or OpenEye."
            )

    def get_rdkit(self, mol: Chem.Mol) -> _T:
        """
        Get the value of this property, given an RDKit molecule.

        Args:
            mol: RDKit molecule.

        Returns:
            Calculated property value.

        Raises:
            NotImplementedError: Property can not be calculated using RDKit.
        """
        if self.composite_func is not None and self.composite_properties is not None:
            return self.composite_func(  # type: ignore
                **{
                    name: property.get_rdkit(mol)
                    for name, property in self.composite_properties.items()
                }
            )
        if self.rdkit_func is None:
            raise NotImplementedError("Property is not calculated using RDKit.")
        return self.rdkit_func(mol)

    def get_chemaxon(self, mol: Molecule) -> _T:
        """
        Get the value of this property, given a Chemaxon molecule.

        Args:
            mol: Chemaxon molecule.

        Returns:
            Calculated property value.

        Raises:
            NotImplementedError: Property can not be calculated using Chemaxon.
        """
        if self.composite_func is not None and self.composite_properties is not None:
            return self.composite_func(  # type: ignore
                **{
                    name: property.get_chemaxon(mol)
                    for name, property in self.composite_properties.items()
                }
            )
        if self.cxmol_func is None:
            raise NotImplementedError("Property is not calculated using Chemaxon.")
        return self.cxmol_func(mol)
