from typing import Literal

from openeye import oechem
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers, rdmolops

from chemutils.rdmol_from_oemol import rdmol_from_oemol


def get_mmff94_forcefield(mol: Chem.RWMol, confId: int = -1) -> AllChem.ForceField:
    """
    Wrapper function around getting an MMFF94 forcefield for a molecule.

    The signature has to match `UFFGetMoleculeForceField`, such that either can be passed into
    ConstrainedEmbed.
    """
    return AllChem.MMFFGetMoleculeForceField(
        mol, AllChem.MMFFGetMoleculeProperties(mol), confId=confId
    )


RDKitForcefield = Literal["UFF", "MMFF94"]

RDKIT_FORCEFIELDS = {
    "UFF": rdForceFieldHelpers.UFFGetMoleculeForceField,
    "MMFF94": get_mmff94_forcefield,
}


def embed_molecule_3d(
    mol: oechem.OEMolBase,
    /,
    *,
    seed: int | None,
    constraint: oechem.OEMolBase | None = None,
    optimize: bool = False,
    forcefield: RDKitForcefield = "UFF",
    timeout: int = 0,
) -> None:
    """
    Embed a molecule in 3D space, generating 3D coordinates using RDKit.

    Optionally, a constraint molecule with existing coordinates can be used.

    Args:
        mol: OpenEye molecule to assign coordinates to.
        constraint: Optional OpenEye molecule to which to constraint the result to. It must be a
            strict subset of the molecule we are embedding.
        seed: If an integer is provided, the random seed to use for embedding. If not provided,
            the result will be non-deterministic.
        optimize: Should forcefield optimization be run after embedding?
        forcefield: Forcefield to use for optimizing and constraining. One of UFF or MMFF94.
        timeout: Timeout in seconds for embedding. If 0, no timeout is used.

    Raises:
        ValueError: Provided constraint molecule does not have 3D coordinates.
        ValueError: Provided constraint molecule is not a subset of the molecule.
        ValueError: If timeout is not 0 (i.e. a timeout is used) and a constraint is provided.
    """
    atom_map: dict[oechem.OEAtomBase, Chem.Atom] = {}

    rdmol = rdmol_from_oemol(mol, atom_map=atom_map)
    rdmol = rdmolops.AddHs(rdmol)

    if constraint is not None:
        if timeout != 0:
            raise ValueError("Timeout is not supported with constrained embedding.")
        _embed_constraint(rdmol, seed=seed, constraint=constraint, forcefield=forcefield)
    else:
        _embed(rdmol, seed=seed, optimize=optimize, forcefield=forcefield, timeout=timeout)

    rdconf = rdmol.GetConformer()

    for oeatom, rdatom in atom_map.items():
        mol.SetCoords(oeatom, list(rdconf.GetAtomPosition(rdatom.GetIdx())))
    mol.SetDimension(3)


def _embed(
    rdmol: Chem.RWMol,
    /,
    *,
    seed: int | None,
    optimize: bool,
    forcefield: RDKitForcefield,
    timeout: int,
) -> None:
    # Can only configure `embedFragmentsSeparately` if we use `EmbedParameters`.
    params = AllChem.EmbedParameters()

    params.randomSeed = seed if seed is not None else -1

    # Make `EmbedParameters` match the defaults of `EmbedMolecule`.
    params.useExpTorsionAnglePrefs = True
    params.useBasicKnowledge = True
    params.useMacrocycleTorsions = True
    params.ETversion = 2
    params.useMacrocycle14config = True
    params.timeout = timeout

    # Use random coords when embedding multiple fragments
    # If not used, the fragments are spaced 100's of angstroms apart
    if len(rdmolops.GetMolFrags(rdmol)) > 1:
        params.embedFragmentsSeparately = False
        params.useRandomCoords = True

    # Infer hybridization, which is required for optimization
    rdmolops.SetHybridization(rdmol)

    rdDistGeom.EmbedMolecule(rdmol, params)

    if optimize:
        if forcefield == "UFF":
            AllChem.UFFOptimizeMolecule(rdmol)
        elif forcefield == "MMFF94":
            AllChem.MMFFOptimizeMolecule(rdmol)
        else:
            raise ValueError("Unknown forcefield.")


def _embed_constraint(
    rdmol: Chem.RWMol,
    /,
    *,
    seed: int | None,
    constraint: oechem.OEMolBase,
    forcefield: RDKitForcefield,
) -> None:
    if constraint.GetDimension() != 3:
        raise ValueError("Constraint molecule does not have 3D coordinates.")
    rdmol_constraint = rdmol_from_oemol(constraint)
    rdmolops.AddHs(rdmol_constraint, False, True)
    AllChem.ConstrainedEmbed(
        rdmol,
        rdmol_constraint,
        randomseed=seed if seed is not None else -1,
        getForceField=RDKIT_FORCEFIELDS[forcefield],
    )
