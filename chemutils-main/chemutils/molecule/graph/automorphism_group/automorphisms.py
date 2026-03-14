from collections import defaultdict

import numpy as np
from openeye import oechem

from .generator import AutomorphismGenerator
from .group import AutomorphismGroup
from .nauty import run_nauty_for_oemol
from .subgroup import AutomorphismSubgroup, as_symmetric_subgroup
from .utils import group_comparator


def get_num_automorphisms(oemol: oechem.OEMolBase, /) -> int:
    """Get the total number of automorphisms for a molecule."""
    group = get_automorphism_group(oemol)
    return group.num_automorphisms


def get_automorphism_group(oemol: oechem.OEMolBase, /) -> AutomorphismGroup:
    nauty = run_nauty_for_oemol(oemol)

    # Determine the fragment (disconnected molecule) each atom belongs to
    # For example, `CCC.CC` will give `[0, 0, 0, 1, 1]`
    _, atom_fragments = oechem.OEDetermineComponents(oemol)
    atom_fragments = np.array(atom_fragments)

    # Group generators by their symmetry classes and fragment
    generators_by_symmetry_classes_and_size = defaultdict(list)
    for generator in nauty.generators:
        multi_fragment = len(set(atom_fragments[generator.indices])) > 1
        symmetry_classes = frozenset(nauty.symmetry_classes[generator.mask])
        generators_by_symmetry_classes_and_size[multi_fragment, symmetry_classes].append(generator)

    # Further split generators into groups that do not overlap
    subgroups = [
        AutomorphismSubgroup(num_atoms=oemol.NumAtoms(), generators=overlapping_generators)
        for generators_by_symmetry_class_and_size in generators_by_symmetry_classes_and_size.values()
        for overlapping_generators in group_comparator(
            items=generators_by_symmetry_class_and_size, comparator=do_generators_overlap
        )
    ]

    # Convert symmetric subgroups into special subclass that is more optimal
    subgroups = [
        symmetric_subgroup
        if (symmetric_subgroup := as_symmetric_subgroup(subgroup)) is not None
        else subgroup
        for subgroup in subgroups
    ]

    return AutomorphismGroup(num_atoms=oemol.NumAtoms(), subgroups=subgroups)


def do_generators_overlap(a: AutomorphismGenerator, b: AutomorphismGenerator) -> bool:
    """Do two automorphisms overlap?"""
    return bool(set(a.indices) & set(b.indices))
