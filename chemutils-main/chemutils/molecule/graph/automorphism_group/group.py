import itertools
import math
from collections.abc import Generator
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .subgroup import AutomorphismSubgroup


@dataclass
class AutomorphismGroup:
    """
    Automorphism group of a molecule.

    The group is not stored explicitly, allowing representation of large automorphism groups
    with trillions of members.

    This allows operations such as randomly sampling automorphisms or calculating RMSDs to be
    performed efficiently.
    """

    num_atoms: int
    """Number of atoms in the molecule."""
    subgroups: list[AutomorphismSubgroup]
    """List of subgroups that make up the full automorphism group."""

    @property
    def num_automorphisms(self) -> int:
        """Get the total number of automorphism subgroups that make up the full automorphism group."""
        result = math.prod(subgroup.num_automorphisms for subgroup in self.subgroups)
        return result

    def iterate_automorphisms(self) -> Generator[npt.NDArray[np.int32], None, None]:
        """
        Iterate over all automorphisms in the group.

        Depending on the group size, this could be very large and slow. Generally, it is better to
        solve your problem without having to explicitly enumerate all automorphisms.
        """
        for subgroup_automorphisms in itertools.product(
            *(subgroup.automorphisms for subgroup in self.subgroups)
        ):
            automorphism = self.identity
            for subgroup_automorphism in subgroup_automorphisms:
                automorphism = automorphism[subgroup_automorphism]
            yield automorphism

    @property
    def num_subgroups(self) -> int:
        """Number of subgroups that make up the full automorphism group."""
        return len(self.subgroups)

    @property
    def identity(self) -> npt.NDArray[np.int32]:
        """Identity automorphism."""
        return np.array(range(self.num_atoms))

    def random_automorphism(self) -> npt.NDArray[np.int32]:
        """Get a uniformly randomly chosen automorphism."""
        automorphism = self.identity
        for subgroup in self.subgroups:
            automorphism = automorphism[subgroup.random_automorphism()]
        return automorphism
