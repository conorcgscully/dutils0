from io import StringIO

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Structure import Structure


def protein_to_pdbstr(protein: Structure) -> str:
    """
    Convert a BioPython structure to its PDB string representation.

    Args:
        protein: BioPython protein structure

    Returns:
        PDB-formatted string representation of the protein
    """
    buffer = StringIO()
    pdb_io = PDBIO()
    pdb_io.set_structure(protein)
    # .save only accepts file as per doc, but any object that has implementation of .write for strings will do
    pdb_io.save(buffer)
    return buffer.getvalue()
