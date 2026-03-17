__all__ = [
    "extract_selected_atoms_as_pdb",
    "gemmi_from_rdkit",
    "gemmi_to_pymol",
    "gemmi_to_rdkit",
    "get_coord_array",
    "mol_coordinates_array",
    "distance",
    "centroid",
    "protein_to_rdkit",
    "pymol_select",
    "read_protein_bytes",
    "read_protein_str",
    "write_protein_bytes",
    "write_protein_str",
    "flag_protein_residues_by_ligand_proximity"
]


def __getattr__(name: str):
    if name == "extract_selected_atoms_as_pdb":
        from .pdb_extract import extract_selected_atoms_as_pdb
        return extract_selected_atoms_as_pdb
    if name == "gemmi_from_rdkit":
        from .convert import gemmi_from_rdkit
        return gemmi_from_rdkit
    if name == "gemmi_to_pymol":
        from .convert import gemmi_to_pymol
        return gemmi_to_pymol
    if name == "gemmi_to_rdkit":
        from .convert import gemmi_to_rdkit
        return gemmi_to_rdkit
    if name == "get_coord_array":
        from .mol3d import get_coord_array
        return get_coord_array
    if name == "distance":
        from .mol3d import distance
        return distance
    if name == "centroid":
        from .mol3d import centroid
        return centroid
    if name == "protein_to_rdkit":
        from .convert import protein_to_rdkit
        return protein_to_rdkit
    if name == "pymol_select":
        from .select import pymol_select
        return pymol_select
    if name == "read_protein_bytes":
        from .read import read_protein_bytes
        return read_protein_bytes
    if name == "read_protein_str":
        from .read import read_protein_str
        return read_protein_str
    if name == "write_protein_bytes":
        from .write import write_protein_bytes
        return write_protein_bytes
    if name == "write_protein_str":
        from .write import write_protein_str
        return write_protein_str
    if name == "flag_protein_residues_by_ligand_proximity":
        from .proximity import flag_protein_residues_by_ligand_proximity
        return flag_protein_residues_by_ligand_proximity
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
