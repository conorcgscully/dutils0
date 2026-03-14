import logging
import pickle
from pathlib import Path

import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.openeye import add_jupyter_oemol_display, initialize_license


def test_no_license(monkeypatch):
    """Test that an EnvironmentError is raised if there is no license file."""
    monkeypatch.setattr(oechem, "OEChemIsLicensed", lambda: False)
    with pytest.raises(EnvironmentError) as e:
        initialize_license()
    assert str(e.value) == "Failed to license OEChem."


def test_add_jupyter_oemol_display():
    """Test that the Jupyter display function is added to `OEMolBase`."""
    add_jupyter_oemol_display()
    assert hasattr(oechem.OEMolBase, "_repr_png_")


def test_error_message_set(caplog):
    with pytest.raises(ValueError):
        oemol_from_smiles("CC(")
    assert [record.msg for record in caplog.records] == [
        "Problem parsing SMILES:",
        "Unclosed branch.",
        "CC(",
        "  ^",
    ]


def test_error_message_debug(caplog):
    caplog.set_level(logging.DEBUG)

    oemol = oemol_from_smiles("[NH2]")

    for atom in oemol.GetAtoms():
        oechem.OECheckAtomValence(atom)

    assert caplog.records[0].msg == "Element: N Charge: +0 Env: 0x000002"
    assert caplog.records[0].levelno == logging.DEBUG


def test_error_message_debug_filter(caplog):
    caplog.set_level(logging.WARNING)

    oemol = oemol_from_smiles("[NH2]")

    for atom in oemol.GetAtoms():
        oechem.OECheckAtomValence(atom)

    assert len(caplog.records) == 0


def test_openeye_serialisation_fix():
    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, "C=C")
    oechem.OEAssignHybridization(oemol)
    assert [atom.GetHyb() for atom in oemol.GetAtoms()] == [2, 2]

    # Old-style serialisation is broken
    with (Path(__file__).parent / "ethylene-serialised-old.pickle").open("rb") as f:
        oemol_old = pickle.load(f)
    assert [atom.GetHyb() for atom in oemol_old.GetAtoms()] == [0, 0]

    # New-style serialisation is correct
    serialised = pickle.dumps(oemol)
    oemol_new = pickle.loads(serialised)
    assert [atom.GetHyb() for atom in oemol_new.GetAtoms()] == [2, 2]
