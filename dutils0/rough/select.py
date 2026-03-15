"""PyMOL selection and object-creation helpers."""

from __future__ import annotations

from typing import Any
# from pymol import cmd


from typing import Any

def pymol_select_to_object(cmd: Any, selection: str, object_name: str) -> str:
    """
    Create a new PyMOL object from a selection by round-tripping through a ChemPy model.
    Returns the name of the new PyMOL object.
    """
    if cmd is None:
        raise TypeError("cmd must be a PyMOL command object")
    if not isinstance(selection, str) or not selection.strip():
        raise ValueError("selection must be a non-empty string")
    if not isinstance(object_name, str) or not object_name.strip():
        raise ValueError("object_name must be a non-empty string")

    get_model = getattr(cmd, "get_model", None)
    load_model = getattr(cmd, "load_model", None)
    count_atoms = getattr(cmd, "count_atoms", None)

    if not callable(get_model) or not callable(load_model) or not callable(count_atoms):
        raise TypeError(
            "cmd must provide callable get_model(), load_model(), and count_atoms()"
        )

    sel = selection.strip()
    name = object_name.strip()

    n = count_atoms(sel)
    if n == 0:
        raise ValueError("selection matches no atoms")

    model = get_model(sel)
    load_model(model, name)

    return name