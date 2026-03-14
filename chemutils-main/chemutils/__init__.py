import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

import warnings

warnings.filterwarnings(
    "ignore", message="The value of the smallest subnormal"
)  # filter out annoying NumPy warnings

try:
    import prody  # type: ignore

    # ProDy has its own terrible implementation of logging...
    prody.LOGGER._logger.setLevel(logging.CRITICAL)
except ImportError:
    pass

from .openeye import (
    add_jupyter_oemol_display,
    fix_openeye_pickling,
    initialize_license,
    setup_custom_error_handling,
)

initialize_license()
add_jupyter_oemol_display()
setup_custom_error_handling()
fix_openeye_pickling()
