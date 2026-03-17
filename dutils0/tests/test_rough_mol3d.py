"""Tests for dutils0.rough.mol3d."""

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dutils0.rough.mol3d import mol_coordinates_array


def _embed_3d(smiles: str, *, seed: int = 0xC0FFEE) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol = Chem.AddHs(mol)
    # Embed one 3D conformer and optionally relax. EmbedMolecule returns 0 on success.
    code = AllChem.EmbedMolecule(mol, randomSeed=seed)
    assert code == 0
    AllChem.UFFOptimizeMolecule(mol)
    assert mol.GetNumConformers() == 1
    assert mol.GetConformer().Is3D()
    return mol


def test_mol_coordinates_array_happy_path_shape_and_dtype():
    mol = _embed_3d("CCO")
    arr = mol_coordinates_array(mol)
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (mol.GetNumAtoms(), 3)
    assert np.issubdtype(arr.dtype, np.floating)
    assert np.isfinite(arr).all()


def test_mol_coordinates_array_values_match_rdkit_positions():
    mol = _embed_3d("CCO")
    expected = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    got = mol_coordinates_array(mol)
    assert got.shape == expected.shape
    assert np.allclose(got, expected)


def test_mol_coordinates_array_selects_explicit_conformer_id():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=2, randomSeed=0xBEEF))
    assert len(ids) == 2
    # EmbedMultipleConfs should produce 3D conformers.
    assert mol.GetConformer(ids[0]).Is3D()
    assert mol.GetConformer(ids[1]).Is3D()

    a0 = mol_coordinates_array(mol, conf_id=ids[0])
    a1 = mol_coordinates_array(mol, conf_id=ids[1])
    assert a0.shape == (mol.GetNumAtoms(), 3)
    assert a1.shape == (mol.GetNumAtoms(), 3)
    # Very likely distinct; but test the invariant: each equals RDKit's own positions.
    assert np.allclose(a0, np.asarray(mol.GetConformer(ids[0]).GetPositions(), dtype=float))
    assert np.allclose(a1, np.asarray(mol.GetConformer(ids[1]).GetPositions(), dtype=float))


def test_mol_coordinates_array_molecule_with_no_conformers_raises():
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    assert mol.GetNumConformers() == 0
    with pytest.raises(ValueError, match="no conformers"):
        mol_coordinates_array(mol)


def test_mol_coordinates_array_invalid_conformer_id_raises():
    mol = _embed_3d("CCO")
    with pytest.raises(ValueError, match="does not exist"):
        mol_coordinates_array(mol, conf_id=999999)

