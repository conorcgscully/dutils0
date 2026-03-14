import numpy as np
import numpy.typing as npt
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder

from .residues import PROTEIN_CODE_ONE_TO_THREE

# fmt: off
OPENFOLD_ATOM_NAMES = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
}
"""Standard ordering of atoms for each residue in the Openfold/Alphafold2 standard."""
# fmt: on


def construct_protein(
    *,
    full_atom_coordinates: npt.NDArray[np.float32],
    sequence: str,
    structure_name: str = "",
    residue_b_factors: npt.NDArray[np.float32] | None = None,
    residue_occupancies: npt.NDArray[np.float32] | None = None,
    residue_numbers: npt.NDArray[np.int32] | None = None,
    residue_insertion_codes: list[str | None] | None = None,
    chain_id: str = "A",
) -> Structure:
    """
    Construct a BioPython Structure object from all atom coordinates and single-letter protein sequences.

    The all atom coordinates are expected to be in standard Openfold/Alphafold2 ordering. The sequence
    is used to assign the positions to each atom:

    ```
    sequence = 'AG...'
    [
        pos1,   # Alanine, resid 1, atom 'N'
        pos2,   # Alanine, resid 1, atom 'CA'
        pos3,   # Alanine, resid 1, atom 'C'
        pos4,   # Alanine, resid 1, atom 'O'
        pos5,   # Alanine, resid 1, atom 'CB'
        pos6,   # Gylcine, resid 2, atom 'N'
        pos7,   # Gylcine, resid 2, atom 'CA'
        pos8,   # Gylcine, resid 2, atom 'C'
        pos9,   # Gylcine, resid 2, atom 'O'
        ...
    ]
    ```

    Args:
        full_atom_coordinates: NumPy array of shape (natoms, 3), where N is the number of atoms.
        sequence: Single-letter protein sequence.
        residue_b_factors: Optional per-residue B-factors for storing additional information.
        residue_occupancies: Optional per-residue occupancies for storing additional information.
        residue_numbers: Optional residue numbers. If not provided, residues will be numbered 1..N.
        residue_insertion_codes: Optional residue insertion codes. If not provided, no insertion codes will be used.
        structure_name: Optional name to attach to the structure.
        chain_id: Optional chain ID to label the structure as.

    Returns:
        BioPython structure that can be written to PDB.

    Raises:
        ValueError: The number of atoms implied by the protein sequence does not match the number
            of coordinates.
    """
    if len(full_atom_coordinates.shape) != 2 or full_atom_coordinates.shape[1] != 3:
        raise ValueError(f"Invalid atom coordinate shape: {full_atom_coordinates.shape}")
    builder = StructureBuilder()
    builder.init_structure(structure_name)
    builder.init_model(0)
    builder.init_chain(chain_id)
    builder.init_seg("")

    atomindex = 0
    for resindex, letter in enumerate(sequence):
        resname = PROTEIN_CODE_ONE_TO_THREE[letter]
        atomnames = OPENFOLD_ATOM_NAMES[resname]
        resid = residue_numbers[resindex] if residue_numbers is not None else resindex + 1
        insertion_code = (
            residue_insertion_codes[resindex] if residue_insertion_codes is not None else None
        )
        builder.init_residue(resname, " ", resid, insertion_code if insertion_code else " ")
        residue_b_factor = residue_b_factors[resindex] if residue_b_factors is not None else 0.0
        residue_occupancy = (
            residue_occupancies[resindex] if residue_occupancies is not None else 1.0
        )
        for atomname in atomnames:
            if atomindex == len(full_atom_coordinates):
                raise ValueError("Mismatch between number of coordinates and protein sequence")
            builder.init_atom(
                name=atomname,
                coord=full_atom_coordinates[atomindex],
                b_factor=residue_b_factor,
                occupancy=residue_occupancy,
                altloc=" ",
                fullname=atomname,
                element=atomname[0],
            )
            atomindex += 1

    if atomindex != len(full_atom_coordinates):
        raise ValueError("Mismatch between number of coordinates and protein sequence")

    return builder.get_structure()
