"""Module for molecular similarity calculation."""

import logging
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Generic, TypeVar

import numpy as np
from openeye import oechem
from rdkit import Chem, RDLogger
from rdkit.Chem import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.rdmol_from_oemol import rdmol_from_oemol

rd_logger = RDLogger.logger()
rd_logger.setLevel(RDLogger.CRITICAL)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class InvalidSmilesError(ValueError):
    """SMILES could not be parsed by the underlying toolkit."""

    pass


def _rdmol_from_smiles(smiles: str) -> Chem.Mol:
    """
    Create an RDKit molecule from a SMILES string, passing through OpenEye toolkit as a bridge.

    Args:
        smiles: SMILES string.

    Returns:
        Chem.Mol: RDKit molecule.
    """
    try:
        return rdmol_from_oemol(oemol_from_smiles(smiles))
    except ValueError as e:
        raise InvalidSmilesError(f"Invalid SMILES string: {smiles}") from e


T_InputMoleculeType = TypeVar("T_InputMoleculeType", str, oechem.OEMolBase)
T_OutputMoleculeType = TypeVar("T_OutputMoleculeType", str, oechem.OEMolBase)


def _fingerprint_from_mol(
    mol: T_InputMoleculeType, /, *, radius: int = 3
) -> DataStructs.ExplicitBitVect | None:
    """
    Generate a Morgan fingerprint from a molecule.

    Args:
        mol: SMILES string or OEChem molecule.
        radius: Radius of the Morgan fingerprint.

    Returns:
        DataStructs.ExplicitBitVect: Morgan fingerprint.
    """
    if isinstance(mol, DataStructs.ExplicitBitVect):
        return mol
    if isinstance(mol, oechem.OEMolBase):
        rdmol = rdmol_from_oemol(mol)
    else:
        try:
            rdmol = _rdmol_from_smiles(mol)
        except InvalidSmilesError:
            logger.warning(f"Failed to generate fingerprint for invalid SMILES: `{mol}`")
            return None
    return GetMorganGenerator(radius=radius).GetFingerprint(rdmol)


def calculate_similarity_in_parallel(
    *,
    query: T_InputMoleculeType,
    targets: Sequence[T_OutputMoleculeType],
    radius: int = 3,
    fingerprints: dict[str, DataStructs.ExplicitBitVect] | None = None,
) -> dict[T_OutputMoleculeType, np.float16 | None]:
    """
    Calculate similarity between a query molecule and a list of target molecules.

    Uses the Tanimoto similarity coefficient and Morgan fingerprints.

    Args:
        query: Query molecule.
        targets: Target molecules.
        radius: Radius of the Morgan fingerprint.
        fingerprints: Precomputed fingerprints for the target molecules.

    Returns:
        Dictionary containing the SMILES string and similarity score.
    """
    query_fp = _fingerprint_from_mol(query, radius=radius)
    if query_fp is None:
        raise InvalidSmilesError(f"Query molecule `{query}` is invalid.")

    if fingerprints is not None:
        target_fps = [fingerprints.get(target) for target in targets]
    else:
        target_fps = [None] * len(targets)

    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(_fingerprint_from_mol, mol, radius=radius): i
            for i, mol in enumerate(targets)
            if target_fps[i] is None
        }

        for future in as_completed(futures):
            idx = futures[future]
            try:
                computed_fp = future.result()
                target_fps[idx] = computed_fp
                if computed_fp is not None and fingerprints is not None:
                    fingerprints[targets[idx]] = computed_fp
            except Exception as e:
                logger.warning(f"Error computing fingerprint for {targets[idx]}: {e}")

    valid_fps: list[DataStructs.ExplicitBitVect] = [fp for fp in target_fps if fp is not None]
    similarities: list[float] = DataStructs.BulkTanimotoSimilarity(query_fp, valid_fps)

    # Keep track of the index of the valid fingerprint to match with the target SMILES
    valid_fingerprint_index = 0
    results: dict[T_OutputMoleculeType, np.float16 | None] = {}

    for i, fp in enumerate(target_fps):
        if fp is None:
            results[targets[i]] = None
        else:
            results[targets[i]] = np.float16(similarities[valid_fingerprint_index])
            valid_fingerprint_index += 1

    return results


@dataclass
class SimilarMolecule(Generic[T_InputMoleculeType]):
    """
    Wraps a molecule and its similarity score.

    Doesn't refer to the query molecule for simplicity and to avoid redundancy across many similar molecules.
    """

    molecule: T_InputMoleculeType
    similarity: np.float16


class InvalidNumResultsError(ValueError):
    """Number of results must be greater than 0."""

    pass


def get_most_similar_molecules(
    *,
    query: T_InputMoleculeType,
    targets: Sequence[T_OutputMoleculeType],
    num_results: int,
    radius: int = 3,
    fingerprints: dict[T_OutputMoleculeType, DataStructs.ExplicitBitVect] | None = None,
) -> list[SimilarMolecule[T_OutputMoleculeType]]:
    """
    Calculate the most similar molecules to a query molecule.

    Args:
        query: Query molecule.
        targets: List of target molecules to compare against.
        num_results: Number of most similar molecules to return.
        radius: Radius of the Morgan fingerprint.
        fingerprints: Precomputed fingerprints for the target molecules.

    Returns:
        List of ranked SimilarMolcule named tuples sorted by similarity and limited to the `num_results`.
    """
    if not targets:
        return []
    if num_results <= 0:
        raise InvalidNumResultsError("num_results must be greater than 0.")
    similarity_dict: dict[T_OutputMoleculeType, np.float16 | None] = (
        calculate_similarity_in_parallel(
            query=query, targets=targets, radius=radius, fingerprints=fingerprints
        )
    )
    return sorted(
        (
            SimilarMolecule(molecule=mol, similarity=similarity)
            for mol, similarity in similarity_dict.items()
            if similarity is not None
        ),
        key=lambda x: x.similarity,
        reverse=True,
    )[:num_results]
