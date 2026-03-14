"""
Compatibility package.

This repo's canonical package is `dutils0`. The `rough` package is provided as a
thin import shim to preserve existing imports like:

    from rough import extract_selected_atoms_as_pdb
"""

from dutils0.rough import extract_selected_atoms_as_pdb

__all__ = ["extract_selected_atoms_as_pdb"]

