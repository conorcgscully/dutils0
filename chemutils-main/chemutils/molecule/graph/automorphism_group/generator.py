from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class AutomorphismGenerator:
    """
    Representation of a generator of the automorphism group.

    Internally, this is represented as a list of disjoint cycles, as this is what is returned by NAUTY.
    """

    num_atoms: int
    """Number of atoms in the molecule."""
    cycles: tuple[tuple[int, ...], ...]
    """List of disjoint cycles that represent the generator."""

    @property
    def order(self) -> int:
        """Order of the generator, which is the number of times it must be applied to return to the identity."""
        return max(len(cycle) for cycle in self.cycles)

    @property
    def automorphism(self) -> npt.NDArray[np.int32]:
        """Numpy array representation of the generator that can be used as an index."""
        automorphism = np.array(range(self.num_atoms), dtype=np.int32)
        for cycle in self.cycles:
            for i in range(len(cycle)):
                automorphism[cycle[i]] = cycle[(i + 1) % len(cycle)]
        return automorphism

    @property
    def indices(self) -> list[int]:
        """List of indices that are affected by the generator."""
        return sorted({index for cycle in self.cycles for index in cycle})

    @property
    def mask(self) -> npt.NDArray[np.bool_]:
        """Boolean mask of the indices that are affected by the generator."""
        mask = np.zeros(self.num_atoms, dtype=bool)
        mask[self.indices] = True
        return mask
