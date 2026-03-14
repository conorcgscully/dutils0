import itertools
from collections.abc import Sequence
from typing import Literal

from openeye import oechem

from chemutils.molecule import (
    create_bond,
    get_atom,
    merge_molecules,
    modify_residue,
    remove_atom,
)
from chemutils.protein._chemical_component_dictionary import get_chemical_component
from chemutils.protein.residues import UnknownCodeHandling, residue_names_from_sequence

SEQUENCE_TYPE_TO_INTERRESIDUE_ATOM_NAMES = {
    "protein": ("C", "N", "OXT"),
    "dna": ("P", "O3'", "OP3"),
    "rna": ("P", "O3'", "OP3"),
}
"""Mapping of sequence type to the start, end and condensed atom names for interresidue bonds."""


def create_polymer(
    *,
    sequence: str | None = None,
    residue_names: Sequence[str] | None = None,
    sequence_type: Literal["protein", "dna", "rna"],
    start_residue_number: int | None = None,
    start_serial_number: int | None = None,
    residue_numbers: Sequence[int] | None = None,
    chain_id: str | None = None,
) -> oechem.OEGraphMol:
    """
    Construct an OpenEye molecule for a protein chain, from either a sequence or list of residue names.

    Args:
        sequence: Single-letter sequence, such as `QYMRTGEGFLLVFAIN`. Can be omitted if `residue_names` is provided.
        residue_names: List of residue names, such as `["GLN", "TYR", "MET", "ARG", ...]`. Can be omitted if `sequence`
            is provided. If both `sequence` and `residue_names`are specified, the two must match.
        sequence_type: Type of polymer to create. Must be one of `"protein"`, `"dna"` or `"rna"`.
        residue_numbers: Array of residue numbers. Must be the same length as the provided sequence. If not provided,
            residues will be numbered sequentially starting from `1`.
        start_residue_number: If `residue_numbers` is not provided, the number to start the sequential numbering of
            residues. Cannot be used when `residue_numbers` is provided.
        start_serial_number: Serial number to start the sequential numbering of atoms. Defaults to `1`.
        chain_id: Chain identifier, such as `"A"`. By default, no chain ID will be assigned.

    Raises:
        ValueError: If neither `sequence` or `residue_names` are provided.
        ValueError: If `sequence` and `residue_names` are not compatible.
        ValueError: If both `start_residue_number` and `residue_numbers` are provided.

    """
    if sequence is not None:
        residue_names_sequence = residue_names_from_sequence(
            sequence,
            sequence_type=sequence_type,
            unknown_code_handling=UnknownCodeHandling.ResidueName,
            nonstandard_code_handling=UnknownCodeHandling.Exception,
            allow_ambiguous=True,
        )
        if residue_names is not None and residue_names != residue_names_sequence:
            raise ValueError("`sequence` and `residue_names` are not compatible.")
        residue_names = residue_names_sequence
    else:
        if residue_names is None:
            raise ValueError("Either `sequence` or `residue_names` must be provided.")

    num_residues = len(residue_names)

    if not residue_numbers:
        start_residue_number = start_residue_number if start_residue_number is not None else 1
        residue_numbers = list(range(start_residue_number, num_residues + start_residue_number))
    elif start_residue_number:
        raise ValueError(
            "Cannot provide `start_residue_number` when `residue_numbers` are provided."
        )

    residues = [get_chemical_component(residue_name) for residue_name in residue_names]

    for residue, residue_number in zip(residues, residue_numbers, strict=True):
        for atom in residue.GetAtoms():
            modify_residue(
                atom,
                residue_number=residue_number,
                chain_id=chain_id,
            )

    atom_mapping: dict[oechem.OEAtomBase, oechem.OEAtomBase] = {}
    bond_mapping: dict[oechem.OEBondBase, oechem.OEBondBase] = {}

    oemol = merge_molecules(*residues, atom_mapping=atom_mapping, bond_mapping=bond_mapping)

    start_atom_name, end_atom_name, condensed_atom_name = SEQUENCE_TYPE_TO_INTERRESIDUE_ATOM_NAMES[
        sequence_type
    ]

    # Remove oxygen atoms that are removed by condensation
    for residue in residues[:-1]:
        atom = get_atom(residue, selection=f'name "{condensed_atom_name}"')
        remove_atom(atom_mapping[atom], bond_handling="hydrogens")

    # Form interresidue bonds
    for residue1, residue2 in itertools.pairwise(residues):
        atom1 = get_atom(residue1, selection=f'name "{start_atom_name}"')
        atom2 = get_atom(residue2, selection=f'name "{end_atom_name}"')
        create_bond(oemol, atom1=atom_mapping[atom1], atom2=atom_mapping[atom2], order=1)

    start_serial_number = start_serial_number if start_serial_number is not None else 1

    # Add atom serial numbers
    for i, atom in enumerate(oemol.GetAtoms()):
        modify_residue(
            atom,
            serial_number=start_serial_number + i,
        )

    return oemol


def create_protein(
    *,
    sequence: str | None = None,
    residue_names: Sequence[str] | None = None,
    start_residue_number: int | None = None,
    start_serial_number: int | None = None,
    residue_numbers: Sequence[int] | None = None,
    chain_id: str | None = None,
) -> oechem.OEGraphMol:
    """
    Construct an OpenEye molecule for a protein chain, from either a sequence or list of residue names.

    Args:
        sequence: Single-letter sequence, such as `QYMRTGEGFLLVFAIN`. Can be omitted if `residue_names` is provided.
        residue_names: List of residue names, such as `["GLN", "TYR", "MET", "ARG", ...]`. Can be omitted if `sequence`
            is provided. If both `sequence` and `residue_names`are specified, the two must match.
        residue_numbers: Array of residue numbers. Must be the same length as the provided sequence. If not provided,
            residues will be numbered sequentially starting from `1`.
        start_residue_number: If `residue_numbers` is not provided, the number to start the sequential numbering of
            residues. Cannot be used when `residue_numbers` is provided.
        start_serial_number: Serial number to start the sequential numbering of atoms. Defaults to `1`.
        chain_id: Chain identifier, such as `"A"`. By default, no chain ID will be assigned.

    Raises:
        ValueError: If neither `sequence` or `residue_names` are provided.
        ValueError: If `sequence` and `residue_names` are not compatible.
        ValueError: If both `start_residue_number` and `residue_numbers` are provided.

    """
    return create_polymer(
        sequence=sequence,
        residue_names=residue_names,
        sequence_type="protein",
        start_residue_number=start_residue_number,
        start_serial_number=start_serial_number,
        residue_numbers=residue_numbers,
        chain_id=chain_id,
    )
