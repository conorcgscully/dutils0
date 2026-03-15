__all__ = ["extract_selected_atoms_as_pdb"]


def __getattr__(name: str):
    if name == "extract_selected_atoms_as_pdb":
        from .pdb_extract import extract_selected_atoms_as_pdb
        return extract_selected_atoms_as_pdb
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
