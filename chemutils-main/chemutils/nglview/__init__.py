"""Code for interfacing with NGLView."""

import importlib
from typing import Any


def load_ipython_extension(ipython: Any) -> None:
    """Register IPython magics."""
    from .ipython import NGLMagics

    has_ngl = importlib.util.find_spec("nglview") is not None

    if not has_ngl:
        raise ImportError("chemutils.nglview requires nglview to be installed.")

    ipython.register_magics(NGLMagics)
