import itertools
import math
import random
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .generator import AutomorphismGenerator


@dataclass
class AutomorphismSubgroup:
    """Subgroup of the automorphism group."""

    num_atoms: int
    """Number of atoms in the molecule."""
    generators: list[AutomorphismGenerator]
    """List of generators that define the subgroup."""

    @property
    def num_automorphisms(self) -> int:
        """Number of automorphisms the subgroup contains (including the identity)."""
        return len(self.automorphisms)

    @property
    def automorphisms(self) -> npt.NDArray[np.int32]:
        """List of all automorphisms in the subgroup."""
        return np.asarray(generate_all_automorphisms_from_generators(self.generators))

    def random_automorphism(self) -> npt.NDArray[np.int32]:
        """Get a uniformly randomly chosen automorphism."""
        return self.automorphisms[random.randint(0, self.num_automorphisms - 1)]  # type: ignore

    @property
    def indices(self) -> list[int]:
        """List of indices that are affected by the subgroup."""
        return sorted(set.union(*(set(generator.indices) for generator in self.generators)))

    @property
    def mask(self) -> npt.NDArray[np.bool_]:
        """Boolean mask of the indices that are affected by the subgroup."""
        mask = np.zeros(self.num_atoms, dtype=bool)
        mask[self.indices] = True
        return mask


@dataclass
class AutomorphismSymmetricSubgroup(AutomorphismSubgroup):
    """
    Representation of a symmetric subgroup with more optimal algorithms.

    Symmetric subgroups are pure permutations of indices or sets of indices, and hence
    the number of automorphisms scales as N!, where N is the degree of the subgroup.

    This class provides more efficient algorithms for both generating all automorphisms,
    and for sampling random automorphisms without calculating N! automorphisms.
    """

    elements: npt.NDArray[np.int32]

    @property
    def degree(self) -> int:
        return len(self.elements)

    @property
    def num_automorphisms(self) -> int:
        return math.factorial(self.degree)

    def random_automorphism(self) -> npt.NDArray[np.int32]:
        indices = self.elements.copy()
        np.random.shuffle(indices)
        automorphism = np.arange(self.num_atoms, dtype=np.int32)
        automorphism[self.elements.flatten()] = indices.flatten()
        return automorphism

    @property
    def automorphisms(self) -> npt.NDArray[np.int32]:
        indices = self.elements.flatten()
        permutations = np.asarray(list(itertools.permutations(range(self.degree))))
        num_permutations = len(permutations)
        automorphisms = np.tile(np.arange(self.num_atoms, dtype=np.int32), (len(permutations), 1))
        new_indices = self.elements[permutations].reshape((num_permutations, len(indices)))
        automorphisms[:, indices] = new_indices
        return automorphisms


def as_symmetric_subgroup(
    subgroup: AutomorphismSubgroup, /
) -> AutomorphismSymmetricSubgroup | None:
    subsets = set()
    for generator in subgroup.generators:
        if generator.order != 2:
            return None
        from_ = []
        to = []
        for cycle in generator.cycles:
            from_.append(cycle[0])
            to.append(cycle[1])
        subsets.add(tuple(from_))
        subsets.add(tuple(to))
    if len(subgroup.generators) == len(subsets) - 1:
        return AutomorphismSymmetricSubgroup(
            num_atoms=subgroup.num_atoms,
            generators=subgroup.generators,
            elements=np.asarray(sorted(subsets)),
        )
    return None


def generate_all_automorphisms_from_generators(
    generators: list[AutomorphismGenerator], /
) -> list[npt.NDArray[np.int32]]:
    """Given a set of automorphism generators, generate all automorphisms (including the identity)."""
    identity = np.array(range(generators[0].num_atoms), dtype=np.int32)
    all_automorphisms = [identity]

    def generate_automorphisms_from_automorphism(automorphism: npt.NDArray[np.int32], /) -> None:
        for generator in generators:
            new_automorphism = automorphism[generator.automorphism]
            if any(
                (new_automorphism == existing_automorphism).all()
                for existing_automorphism in all_automorphisms
            ):
                continue
            all_automorphisms.append(new_automorphism)
            generate_automorphisms_from_automorphism(new_automorphism)

    generate_automorphisms_from_automorphism(identity)

    return all_automorphisms
