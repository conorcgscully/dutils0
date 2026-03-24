"""Convert Gemmi Structure or Model to RDKit ROMol or PyMOL object via PDB."""

from __future__ import annotations

from typing import Any

import gemmi
from rdkit import Chem

from .read import read_protein_str
from .write import write_protein_str

_TEMP_PYMOL_OBJECT = "_dutils0_protein_to_rdkit"


def protein_to_rdkit(
    model_or_structure: gemmi.Model | gemmi.Structure,
    *,
    cmd: Any = None,
) -> "rdkit.Chem.rdchem.Mol":
    """
    Convert a Gemmi Model or Structure to an RDKit ROMol via PDB with PyMOL normalization.

    Pipeline: Model → PDB string (write_protein_str) → PyMOL read_pdbstr →
    PyMOL get_pdbstr → RDKit MolFromPDBBlock. PyMOL normalizes the PDB so RDKit
    receives well-formed input. No file I/O; in-memory only.

    Args:
        model_or_structure: A Gemmi Model to convert, or a Structure (first model is used).
        cmd: Optional PyMOL command object. If None, uses global pymol.cmd.

    Returns:
        RDKit molecule (ROMol). Unsanitized; hydrogens kept. Callers may call
        SanitizeMol if needed.

    Raises:
        TypeError: model_or_structure is not a gemmi.Model or gemmi.Structure.
        ValueError: Structure has no models; or PyMOL/RDKit failed.
            Exceptions from write_protein_str are propagated unchanged.
    """
    if isinstance(model_or_structure, gemmi.Structure):
        if len(model_or_structure) < 1:
            raise ValueError("Structure has no models")
        model = model_or_structure[0].clone()
    elif isinstance(model_or_structure, gemmi.Model):
        model = model_or_structure
    else:
        raise TypeError(
            f"model_or_structure must be a gemmi.Model or gemmi.Structure, got {type(model_or_structure).__name__}"
        )
    # Build a single-model Structure for write_protein_str (clone to avoid mutating caller's model).
    st = gemmi.Structure()
    st.add_model(model.clone())
    pdb_str = write_protein_str(st, "pdb")
    if cmd is None:
        import pymol
        cmd = pymol.cmd
    try:
        try:
            cmd.read_pdbstr(pdb_str, _TEMP_PYMOL_OBJECT)
        except Exception as e:
            raise ValueError(f"PyMOL failed to load PDB: {e}") from e
        try:
            pymol_pdb_str = cmd.get_pdbstr(_TEMP_PYMOL_OBJECT)
        except Exception as e:
            raise ValueError(f"PyMOL failed to write PDB string: {e}") from e
        mol = Chem.MolFromPDBBlock(
            pymol_pdb_str, sanitize=False, removeHs=False
        )
        if mol is None:
            raise ValueError(
                "RDKit failed to parse PDB block (MolFromPDBBlock returned None)"
            )
        return mol
    finally:
        _delete_temp_pymol_object(cmd)


def _delete_temp_pymol_object(cmd: Any) -> None:
    """Delete the temporary PyMOL object; swallow errors so they don't mask the original."""
    try:
        cmd.delete(_TEMP_PYMOL_OBJECT)
    except Exception:
        pass


def gemmi_to_rdkit(structure: gemmi.Structure) -> "rdkit.Chem.rdchem.Mol":
    """
    Convert a Gemmi Structure to an RDKit ROMol via PDB-mediated conversion (Gemmi → PDB string → RDKit).

    Uses write_protein_str for serialization; RDKit infers the molecule from the PDB text. The returned
    mol may be unsanitized (sanitize=False); callers that need sanitization can use RDKit's SanitizeMol.

    Args:
        structure: The Gemmi Structure to convert (model → chain → residue → atom).

    Returns:
        RDKit molecule (ROMol). May be unsanitized; suitable for structural use (e.g. rendering, coordinates).

    Raises:
        TypeError: structure is not a gemmi.Structure.
        ValueError: RDKit failed to parse the PDB block (MolFromPDBBlock returned None). Exceptions from
            write_protein_str (e.g. serialization failure) are propagated unchanged.

    Example:
        mol = gemmi_to_rdkit(st)
    """
    if not isinstance(structure, gemmi.Structure):
        raise TypeError(
            f"structure must be a gemmi.Structure, got {type(structure).__name__}"
        )
    pdb_str = write_protein_str(structure, "pdb")
    mol = Chem.MolFromPDBBlock(pdb_str, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(
            "RDKit failed to parse PDB block (MolFromPDBBlock returned None)"
        )
    return mol


def gemmi_from_rdkit(mol: "rdkit.Chem.rdchem.Mol") -> gemmi.Structure:
    """
    Convert an RDKit ROMol to a Gemmi Structure via PDB-mediated conversion (ROMol → PDB string → Gemmi).

    Uses Chem.MolToPDBBlock for serialization; read_protein_str parses the PDB text. The ROMol should
    have a conformer (3D coordinates) for meaningful output; atom naming, residue info, and chain IDs
    in the result may differ from the original if the mol was produced by gemmi_to_rdkit.

    Args:
        mol: The RDKit molecule (ROMol) to convert.

    Returns:
        gemmi.Structure with at least one model (model → chain → residue → atom).

    Raises:
        TypeError: mol is not an RDKit Mol (ROMol).
        ValueError: MolToPDBBlock returned None or empty string. Exceptions from read_protein_str
            (e.g. parse failure, no model) are propagated unchanged.

    Example:
        st = gemmi_from_rdkit(mol)
    """
    if mol is None or not isinstance(mol, Chem.Mol):
        raise TypeError(
            f"mol must be an RDKit Mol (ROMol), got {type(mol).__name__}"
        )
    pdb_str = Chem.MolToPDBBlock(mol)
    if pdb_str is None or not pdb_str.strip():
        raise ValueError(
            "RDKit MolToPDBBlock returned None or empty string"
        )
    return read_protein_str(pdb_str)


def gemmi_to_pymol(
    structure: gemmi.Structure,
    object_name: str,
    *,
    cmd: Any = None,
) -> Any:
    """
    Convert a Gemmi Structure to a PyMOL-loaded object via PDB string (Gemmi → PDB → PyMOL).

    Uses write_protein_str for serialization; pymol.cmd.read_pdbstr loads the PDB text into PyMOL.
    When PyMOL is not installed and cmd is not provided, ImportError may be raised at call time.

    Args:
        structure: The Gemmi Structure to convert (model → chain → residue → atom).
        object_name: Name for the created/updated object in PyMOL. Must be non-empty after stripping.
        cmd: Optional PyMOL command object. If None, uses global pymol.cmd (requires PyMOL installed).

    Returns:
        The PyMOL cmd object that was used (the one passed in or pymol.cmd).

    Raises:
        TypeError: structure is not a gemmi.Structure.
        ValueError: object_name is not a non-empty string, or PyMOL failed to load the PDB.
            Exceptions from write_protein_str (e.g. serialization failure) are propagated unchanged.

    Example:
        cmd = gemmi_to_pymol(st, "mymol")
    """
    if not isinstance(structure, gemmi.Structure):
        raise TypeError(
            f"structure must be a gemmi.Structure, got {type(structure).__name__}"
        )
    if not isinstance(object_name, str) or not object_name.strip():
        raise ValueError("object_name must be a non-empty string")
    name = object_name.strip()
    pdb_str = write_protein_str(structure, "pdb")
    if cmd is None:
        import pymol

        cmd = pymol.cmd
    try:
        cmd.read_pdbstr(pdb_str, name)
    except Exception as e:
        raise ValueError(f"PyMOL failed to load PDB: {e}") from e
    return cmd
