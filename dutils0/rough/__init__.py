__all__ = [
    "extract_selected_atoms_as_pdb",
    "gemmi_from_rdkit",
    "gemmi_to_pymol",
    "gemmi_to_rdkit",
    "protein_to_rdkit",
    "pymol_select",
    "read_protein_bytes",
    "read_protein_str",
    "write_protein_bytes",
    "write_protein_str",
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
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
