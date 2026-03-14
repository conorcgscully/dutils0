"""
Code for running NAUTY, a library for computing generators of automorphism groups.

The existing pynauty python wrapper is GPL, so we have to write our own wrapper around
the NAUTY library (which is available under an Apache 2.0 license through conda).

Documentation for the NAUTY tool can be found at:
    https://users.cecs.anu.edu.au/~bdm/nauty/nug28.pdf

We have to wrap the `dreadnaut` CLI tool, which is the easiest way to interact with the
NAUTY library without generating bindings for the C library.
"""

import re
import subprocess
import tempfile
from dataclasses import dataclass
from shutil import which

import fsutils as fs
import numpy as np
import numpy.typing as npt
from openeye import oechem

from .generator import AutomorphismGenerator
from .utils import group_comparator

OVERVIEW_REGEX = re.compile(
    r"^(?P<orbits>\d+) orbits?; grpsize=(?P<num>\d+\.?\d*e?\d*); \d+ gens?; \d+ nodes?; maxlev=\d+$"
)


@dataclass
class NautyOutput:
    """Wrapper around the output of `pynauty.autgrp`, which calculates Automorphisms."""

    generators: list[AutomorphismGenerator]
    """List of generators of the automorphism group."""
    num_automorphisms: int | None
    """Total number of automorphisms in the group, or `None` if it was too large."""
    symmetry_classes: npt.NDArray[np.int32]
    """Array of length `num_vertices` where each element is the index of the symmetry class of the corresponding vertex."""
    num_symmetry_classes: int
    """Number of unique symmetry classes."""


def run_nauty_for_oemol(oemol: oechem.OEGraphMol, /) -> NautyOutput:
    """Get the automorphisms of a graph using PyNAUTY."""
    input = generate_dreadnaut_input_from_oemol(oemol)
    output = run_dreadnaut(input=input)

    generator_cycles: list[tuple[tuple[int, ...], ...]] = []

    for line in output.splitlines():
        if line.startswith("("):
            generator_cycles.append(parse_cycles(line))

    overview = [
        match for line in output.splitlines() if (match := OVERVIEW_REGEX.match(line)) is not None
    ]
    if len(overview) == 0:
        raise ValueError("Could not find overview line in dreadnaut output")

    num_orbits = int(overview[0].group("orbits"))
    num_automorphisms = (
        int(overview[0].group("num")) if "." not in overview[0].group("num") else None
    )

    orbits = np.array(range(oemol.NumAtoms()))

    # Calculate the orbits (symmetry classes)
    # Two vertices are in the same orbit if they can be interchanged by an automorphism
    for cycles in generator_cycles:
        for disjoint_cycle in cycles:
            idx1 = disjoint_cycle[0]
            for idx2 in disjoint_cycle[1:]:
                if orbits[idx1] != orbits[idx2]:
                    color_to = orbits[idx1]
                    color_from = orbits[idx2]
                    for i in range(len(orbits)):
                        if orbits[i] == color_from:
                            orbits[i] = color_to

    assert len(np.unique(orbits)) == num_orbits

    return NautyOutput(
        generators=[
            AutomorphismGenerator(num_atoms=oemol.NumAtoms(), cycles=cycles)
            for cycles in generator_cycles
        ],
        num_automorphisms=num_automorphisms,
        num_symmetry_classes=num_orbits,
        symmetry_classes=orbits,
    )


def generate_dreadnaut_input_from_oemol(oemol: oechem.OEMolBase, /) -> str:
    """
    Generate the input script that is provided to the `dreadnaut` CLI tool.

    See Section 2 of the following document for more information:
        https://users.cecs.anu.edu.au/~bdm/nauty/nug28.pdf
    """
    lines = []
    # Use sparse mode for nauty
    lines.append("As")
    # No line limit. Prevents line wrapping that prevent reading the orbits correctly
    lines.append("l=0")
    # Set number of vertices
    lines.append(f"n={oemol.NumAtoms()} g")
    # Generate the adjacency list
    for atom in oemol.GetAtoms():
        lines.append(
            f" {atom.GetIdx()} : {' '.join(str(neigh.GetIdx()) for neigh in atom.GetAtoms() if neigh.GetIdx() > atom.GetIdx())} ;"
        )

    # Partition the atoms by their atomic number and aromaticity
    # This tells NAUTY which vertices can be interchanged
    atom_colour = [(atom.GetAtomicNum(), atom.IsAromatic()) for atom in oemol.GetAtoms()]
    partitions = group_comparator(
        items=list(range(oemol.NumAtoms())),
        comparator=lambda a, b: atom_colour[a] == atom_colour[b],
    )
    lines.append(
        "f=" + "[" + "|".join(",".join(str(i) for i in group) for group in partitions) + "]"
    )

    # Run NAUTY and print out automorphisms
    lines.append("x")
    return "\n".join(lines)


def run_dreadnaut(*, input: str) -> str:
    if which("dreadnaut") is None:
        raise RuntimeError(
            "`dreadnaut` CLI not found. Have you installed NAUTY using `conda install -c conda-forge nauty`?"
        )
    with tempfile.TemporaryDirectory() as tempdir:
        input_path = f"{tempdir}/input.txt"
        fs.write_str(input_path, contents=input)
        # call dreadnaut
        output = subprocess.check_output(["dreadnaut"], input=f"< {input_path}", text=True)
        return output


def parse_cycles(line: str, /) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(v) for v in orbit_str.split()) for orbit_str in line[1:-1].split(")("))
