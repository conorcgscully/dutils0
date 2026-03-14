from collections.abc import Collection

import numpy as np
from openeye import oechem

from chemutils.maths.convert import as_coordinates
from chemutils.maths.types import Vectors
from chemutils.protein.mapping import get_alpha_carbon_mapping_by_residue_number

GDT_TS_CUTOFFS = (1.0, 2.0, 4.0, 8.0)


def get_unaligned_gdt_ts(src: oechem.OEMol, dest: oechem.OEMol, /) -> float:
    """
    Computes GDT_TS as defined in https://en.wikipedia.org/wiki/Global_distance_test .

    This implementation does not align the two proteins before computing the GDT_TS score,
    which means it does not return the best possible GDT_TS score.
    """
    alpha_carbon_mapping = get_alpha_carbon_mapping_by_residue_number(src, dest)

    src_ca_coords = np.array([as_coordinates(atom) for atom in alpha_carbon_mapping])
    dest_ca_coords = np.array([as_coordinates(atom) for atom in alpha_carbon_mapping.values()])

    return _unaligned_gdt(src_ca_coords, dest_ca_coords, cutoffs=GDT_TS_CUTOFFS)


def _unaligned_gdt(
    src_ca_coords: Vectors, dest_ca_coords: Vectors, /, cutoffs: Collection[float]
) -> float:
    """Global distance test (GDT) for two sets of coordinates, without alignment."""
    if src_ca_coords.shape != dest_ca_coords.shape:
        raise ValueError(f"Shape mismatch: {src_ca_coords.shape} != {dest_ca_coords.shape}")

    distances = np.linalg.norm(src_ca_coords - dest_ca_coords, ord=2, axis=-1)

    gdt_results = []
    for cutoff in cutoffs:
        gdt_results.append(np.mean(distances < cutoff))

    return np.mean(gdt_results)  # type: ignore
