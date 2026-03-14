#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import annotations

import math
import os
from functools import lru_cache

import fsutils as fs
from openeye import oechem
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemutils.props.property import MolecularProperty
from chemutils.rdmol_from_oemol import as_rdmol, clean_rdmol_hydrogens


@lru_cache
def get_fragment_scores() -> dict[int, float]:
    data = fs.read_pickle(f"{os.path.dirname(__file__)}/fpscores.pkl.gz")
    fragment_scores = {}
    for i in data:
        for j in range(1, len(i)):
            fragment_scores[i[j]] = float(i[0])
    return fragment_scores


def get_synthetic_accessibility_score(mol: Chem.Mol | oechem.OEMolBase) -> float:
    """
    Get the synthetic accessibility score of a molecule, based on Ertl and Schuffenhauer (2009).

    This estimates the ease of synthesis of a molecule.

    This implementation is based on https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py.

    Args:
        mol: RDKit molecule.

    Returns:
        Synthetic Accessibility Score for the molecule.
    """
    mol = as_rdmol(mol)
    mol = clean_rdmol_hydrogens(mol)

    fscores = get_fragment_scores()

    # fragment score
    fingerprint = rdMolDescriptors.GetMorganFingerprint(mol, radius=2)
    nonzero_fingerprint = fingerprint.GetNonzeroElements()
    score1 = sum(fscores.get(bitId, -4) * v for bitId, v in nonzero_fingerprint.items()) / sum(
        nonzero_fingerprint.values()
    )

    # features score
    num_heavy_atoms = sum(atom.GetSymbol() != "H" for atom in mol.GetAtoms())
    num_chiral_centers = len(
        Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    )
    num_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    num_bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    num_macrocycles = sum(len(ring) > 8 for ring in mol.GetRingInfo().AtomRings())

    penalty_size = num_heavy_atoms**1.005 - num_heavy_atoms
    penalty_stereo = math.log10(num_chiral_centers + 1)
    penalty_spiro = math.log10(num_spiro + 1)
    penalty_bridgehead = math.log10(num_bridgehead + 1)

    # This differs from the paper, which defines it as math.log10(num_macrocycles+1)
    # This generates better results when 2 or more macrocycles are present
    penalty_macrocycle = math.log10(2) if num_macrocycles else 0

    score2 = (
        -penalty_size - penalty_stereo - penalty_spiro - penalty_bridgehead - penalty_macrocycle
    )

    # correction for the fingerprint density
    # not in the original publication, added in version 1.1
    # to make highly symmetrical molecules easier to synthetise
    score3 = 0.0
    if num_heavy_atoms > len(nonzero_fingerprint):
        score3 = math.log(float(num_heavy_atoms) / len(nonzero_fingerprint)) * 0.5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    min = -4.0
    max = 2.5
    sascore = 11.0 - (sascore - min + 1) / (max - min) * 9.0
    # smooth the 10-end
    if sascore > 8.0:
        sascore = 8.0 + math.log(sascore + 1.0 - 9.0)
    if sascore > 10.0:
        sascore = 10.0
    elif sascore < 1.0:
        sascore = 1.0

    return sascore  # type: ignore


SyntheticAccessibilityScore = MolecularProperty(
    id="synthetic_accessibility", rdkit_func=get_synthetic_accessibility_score
)
