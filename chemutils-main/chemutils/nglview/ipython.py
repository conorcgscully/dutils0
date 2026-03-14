"""IPython magick for using OpenMM."""

from typing import Any

import nglview
from IPython.core.magic import Magics, line_magic, magics_class

from .widget import get_nglview_structure


@magics_class
class NGLMagics(Magics):
    """IPython magic for use with NGL."""

    @line_magic  # type: ignore[misc]
    def nglview(self, line: str) -> Any:
        """Display the given object using nglview."""
        return nglview.NGLWidget(get_nglview_structure(self.shell.ev(line)))  # type: ignore
