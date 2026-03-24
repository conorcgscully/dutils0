"""PyMOL selection and object-creation helpers."""

from __future__ import annotations

from typing import Any
# from pymol import cmd


from rdkit import Chem


def select_residue_as_mol(mol: Chem.Mol, chain_resi: str) -> Chem.Mol:
    """
    Return a new RDKit molecule containing only the selected residue.

    The selector must be strictly in the form ``CHAIN:RESI`` (exactly one ``:``),
    where CHAIN is non-empty and RESI is an integer.

    Design notes:
    - Query string: `chain <CHAIN> and resi <RESI>` is the most direct rdsl expression.
    - Assumes the RDKit molecule is PDB-derived with `AtomPDBResidueInfo` present.
    - Caveats: insertion codes / altlocs are not expressible via `CHAIN:RESI`.
    """
    if not isinstance(chain_resi, str):
        raise TypeError(f"chain_resi must be a str like 'A:42', got {type(chain_resi).__name__}")

    s = chain_resi.strip()
    if s.count(":") != 1:
        raise ValueError(f"selector must be in form CHAIN:RESI (exactly one ':'), got {chain_resi!r}")

    chain_raw, resi_raw = s.split(":", 1)
    chain = chain_raw.strip()
    if not chain:
        raise ValueError(f"selector must have a non-empty chain id before ':', got {chain_resi!r}")

    resi_str = resi_raw.strip()
    try:
        resi = int(resi_str)
    except ValueError as e:
        raise ValueError(f"selector residue must be an integer, got {resi_str!r}") from e

    try:
        from rdsl import select_molecule
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "rdsl is required; ensure vendored `rdsl/src` is importable or install rdsl."
        ) from e

    expr = f"chain {chain} and resi {resi}"
    result = select_molecule(mol, expr)

    submol = getattr(result, "mol", None)
    if submol is None:
        raise ValueError(f"no atoms match selector {chain_resi!r} (rdsl expr: {expr!r})")
    if not hasattr(submol, "GetNumAtoms"):
        raise TypeError(f"rdsl returned unexpected result.mol type: {type(submol).__name__}")
    if submol.GetNumAtoms() == 0:
        raise ValueError(f"rdsl returned an empty molecule for selector {chain_resi!r} (expr: {expr!r})")

    return Chem.Mol(submol)


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