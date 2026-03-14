import logging
import os
import traceback
from typing import Any

import fsutils as fs
from openeye import oechem


def initialize_license() -> None:
    """Check if OpenEye is licensed, and download a license from S3 if not."""
    if (
        not (
            "OE_LICENSE" in os.environ or "OE_DIR" in os.environ or os.path.exists("oe_license.txt")
        )
        or not oechem.OEChemIsLicensed()  # don't even bother calling oechem if there's no licence
    ):  # pragma: no cover
        try:
            license = fs.read_bytes("s3://charmtx-datalake/licenses/oe_license.txt").decode()
            oechem.OEAddLicenseData(license)

        except Exception as e:
            raise OSError("Failed to download OEChem license from S3.") from e
        if not oechem.OEChemIsLicensed():
            raise OSError("Failed to license OEChem.")


def add_jupyter_oemol_display() -> None:
    """Add a PNG representation of `OEMolBase` to Jupyter, based on converting to RDKit first."""
    try:
        from chemutils.draw import (
            DEFAULT_HIGHLIGHT,
            JUPYTER_HEIGHT,
            JUPYTER_WIDTH,
            Highlight,
            draw_molecule,
        )
        from chemutils.draw.color import get_colors
    except ImportError:
        return

    def default_draw(
        mol: oechem.OEMolBase,
        highlights: Highlight | list[Highlight] | None = None,
    ) -> bytes:
        return draw_molecule(
            mol, width=JUPYTER_WIDTH, height=JUPYTER_HEIGHT, format="png", highlight=highlights
        )

    def _repr_png_oemol_(self: oechem.OEMolBase) -> bytes:
        return default_draw(self)

    oechem.OEMolBase._repr_png_ = _repr_png_oemol_

    def _repr_png_oeatom_(self: oechem.OEAtomBase) -> bytes:
        return default_draw(
            self.GetParent(), highlights=Highlight(atoms=self, color=DEFAULT_HIGHLIGHT)
        )

    oechem.OEAtomBase._repr_png_ = _repr_png_oeatom_

    def _repr_png_oeatomiter_(self: oechem.OEAtomIter) -> bytes | None:
        atoms = list(self)
        if len(atoms) == 0:
            return None
        return default_draw(
            atoms[0].GetParent(), highlights=Highlight(atoms=atoms, color=DEFAULT_HIGHLIGHT)
        )

    oechem.OEAtomIter._repr_png_ = _repr_png_oeatomiter_

    def _repr_png_oebond_(self: oechem.OEBondBase) -> bytes:
        return default_draw(
            self.GetParent(), highlights=Highlight(bonds=self, color=DEFAULT_HIGHLIGHT)
        )

    oechem.OEBondBase._repr_png_ = _repr_png_oebond_

    def _repr_png_oebonditer_(self: oechem.OEAtomIter) -> bytes | None:
        bonds = list(self)
        if len(bonds) == 0:
            return None
        return default_draw(
            bonds[0].GetParent(), highlights=Highlight(bonds=bonds, color=DEFAULT_HIGHLIGHT)
        )

    oechem.OEBondIter._repr_png_ = _repr_png_oebonditer_

    def _repr_png_oeatombondset_(self: oechem.OEAtomBondSet) -> bytes | None:
        if self.NumAtoms() == 0 and self.NumBonds() == 0:
            return None
        return default_draw(
            self.GetParent(), highlights=Highlight(atoms=self, bonds=self, color=DEFAULT_HIGHLIGHT)
        )

    oechem.OEAtomBondSet._repr_png_ = _repr_png_oeatombondset_

    def _repr_png_oeatombondsetvector_(self: oechem.OEAtomBondSetVector) -> bytes | None:
        if len(self) == 0:
            return None
        colors = get_colors(len(self))
        return default_draw(
            self[0].GetParent(),
            highlights=[
                Highlight(atoms=atombondset, bonds=atombondset, color=colors[i])
                for i, atombondset in enumerate(self)
            ],
        )

    oechem.OEAtomBondSetVector._repr_png_ = _repr_png_oeatombondsetvector_


def fix_openeye_pickling() -> None:
    """Change `OEGraphMol` pickling to include properties not serialized to OEB."""

    def __getstate__oemol(self: oechem.OEGraphMol) -> dict[str, Any]:
        # Override OpenEye's default serialisation to bytes
        # Serialise a dictionary with the base OpenEye serialisation and additional data
        # that is not serialised as additional keys
        oeb_bytes = oechem.OEWriteMolToBytes(".oeb", self)
        atom_chiral = [atom.IsChiral() for atom in self.GetAtoms()]
        bond_chiral = [bond.IsChiral() for bond in self.GetBonds()]
        hybridisation = [atom.GetHyb() for atom in self.GetAtoms()]
        return {
            "bytes": oeb_bytes,
            "atom_chiral": atom_chiral,
            "bond_chiral": bond_chiral,
            "atom_hyb": hybridisation,
        }

    def __setstate__oemol(self: oechem.OEGraphMol, state: Any) -> None:
        mol = oechem.OEGraphMol()
        # Backwards compatibility with loading pure OpenEye state, which is just a bytes object.
        if isinstance(state, bytes):
            oechem.OEReadMolFromBytes(mol, ".oeb", state)
            self.__init__(mol)
        else:
            # Read our custom serialisation state.
            oechem.OEReadMolFromBytes(mol, ".oeb", state["bytes"])
            self.__init__(mol)
            for atom, chirality in zip(self.GetAtoms(), state["atom_chiral"], strict=True):
                atom.SetChiral(chirality)
            for atom, hybridisation in zip(self.GetAtoms(), state["atom_hyb"], strict=True):
                atom.SetHyb(hybridisation)
            for bond, chirality in zip(self.GetBonds(), state["bond_chiral"], strict=True):
                bond.SetChiral(chirality)
        del mol

    oechem.OEGraphMol.__getstate__ = __getstate__oemol
    oechem.OEGraphMol.__setstate__ = __setstate__oemol


OE_ERROR_LEVEL_TO_LOGGING = {
    oechem.OEErrorLevel_Debug: logging.DEBUG,
    oechem.OEErrorLevel_Info: logging.INFO,
    oechem.OEErrorLevel_Warning: logging.WARNING,
    oechem.OEErrorLevel_Error: logging.ERROR,
    oechem.OEErrorLevel_Fatal: logging.FATAL,
}


class CustomErrorHandling(oechem.OEErrorHandlerImplBase):  # type: ignore[misc]
    """Custom OpenEye error handler that writes to logging."""

    def Msg(self, level: int, arg0: str) -> None:
        # Get the preceeding frame to print the warning as if it were thrown on that line
        frame_summary = traceback.extract_stack(limit=2)[0]
        function_name = frame_summary.name
        line_no = frame_summary.lineno
        module = frame_summary.filename.split("openeye")[-1].removeprefix("/").removesuffix(".py")
        openeye_logger = logging.getLogger(f"openeye.{module}.{function_name}")
        level = OE_ERROR_LEVEL_TO_LOGGING.get(level, logging.INFO)

        if not openeye_logger.isEnabledFor(level):
            return

        record = openeye_logger.makeRecord(
            name=f"openeye.{module}.{function_name}",
            level=level,
            exc_info=None,
            args={},
            fn=frame_summary.filename,
            func=function_name,
            lno=line_no if line_no else 0,
            msg=arg0.rstrip("\n"),
        )
        openeye_logger.handle(record)


# Has to be in global scope to prevent segfaults from it being deleted

OPENEYE_ERROR_HANDLER = CustomErrorHandling()


def setup_custom_error_handling() -> None:
    """Setup the error handling for OpenEye."""
    oechem.OEThrow.SetHandlerImpl(OPENEYE_ERROR_HANDLER, owned=False)
