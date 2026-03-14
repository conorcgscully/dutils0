import pytest
from openeye import oechem, oeomega

from chemutils.molecule import (
    assert_mol_equal,
    get_coordinates,
    standardise_oemol,
    std_oemol_from_smiles,
)
from chemutils.molecule.conformer import (
    OEOmegaError,
    generate_conformers_omega,
    generate_single_conformer,
)
from chemutils.molecule.equal import AllAtomProperties, AtomProperty


def test_generate_conformers_omega_multiple():
    opts = oeomega.OEOmegaOptions()
    opts.SetMaxConfs(5)

    oemol = std_oemol_from_smiles("CCCCCCCCCCCC")

    confs = generate_conformers_omega(oemol, options=opts)
    standardise_oemol(confs)

    assert isinstance(confs, oechem.OEMCMolBase)
    assert confs.NumConfs() == 5
    assert get_coordinates(confs).shape == (5, confs.NumAtoms(), 3)
    assert get_coordinates(confs.GetConf(oechem.OEHasConfIdx(0))).shape == (confs.NumAtoms(), 3)
    assert_mol_equal(oemol, confs, atom_properties=AllAtomProperties & ~AtomProperty.Coordinates)


def test_generate_conformers_omega_error():
    oemol = std_oemol_from_smiles("C(F)(Cl)(Br)I")

    with pytest.raises(OEOmegaError) as e:
        _ = generate_conformers_omega(oemol)

    assert str(e.value) == "UnspecifiedStereo"


def test_generate_single_conformer():
    oemol = std_oemol_from_smiles("CCC")

    confs = generate_single_conformer(oemol)
    standardise_oemol(confs)

    assert isinstance(confs, oechem.OEConfBase)
    assert get_coordinates(confs).shape == (confs.NumAtoms(), 3)
    assert_mol_equal(oemol, confs, atom_properties=AllAtomProperties & ~AtomProperty.Coordinates)
