![workflow](https://github.com/CHARM-Tx/chemutils/actions/workflows/full_run.yml/badge.svg)
[![codecov](https://codecov.io/github/CHARM-Tx/chemutils/branch/main/graph/badge.svg?token=IN9FIRIC0X)](https://codecov.io/github/CHARM-Tx/chemutils)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# chemutils

`chemutils` is a collection of utility methods and common functionality relating to chemistry. This repo is for short functions with well-defined behaviour - more complex functions probably belong in the repo they are required.

## Installing

To install the latest version of chemutils, run:

`pip install ch-chemutils`

## Environment Set Up

To set up a development environment for chemutils using conda:

```bash
conda create -n chemutils python=3.10
conda activate chemutils
pip install -e .[tests]
```

## Development

This repo uses `ruff` for linting and formatting, and `mypy` for type checking. To run these, use:

```
ruff . --fix
ruff format .
mypy
```

## Tests

Tests can be run using

```bash
python -m pytest
```

# Molecular Properties

`chemutils` provides the following molecular properties:

| Property in `chemutils.props`        | ID                                  | Description                                                                                            |
| ------------------------------------ | ----------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `MolecularWeight`                    | `'molecular_weight'`                | Molecular weight, according to OpenEye                                                                 |
| `NumberHeavyAtoms`                   | `'num_heavy_atoms'`                 | Number of non-hydrogen atoms                                                                           |
| `NumberChiralAtoms`                  | `'num_chiral_atoms'`                | Number of chiral atoms (tetrahedral stereocenters)                                                     |
| `SpacialScore`                       | `'spacial_score'`                   | Spacial Score (measure of molecular complexity)                                                        |
| `NumberSpiroAtoms`                   | `'num_spiro_atoms'`                 | Number of spiro-center atoms (tetrahedral stereocenters)                                               |
| `NumberRotatableBonds`               | `'num_rotatable_bonds'`             | Number of rotatable bonds as defined by OpenEye                                                        |
| `NumberAromaticRings`                | `'num_aromatic_rings'`              | Number of aromatic rings as defined by OpenEye                                                         |
| `QuantitativeEstimationDrugLikeness` | `'qed'`                             | Quantitative Estimation of Drug-Likeness, as calculated by RDKit                                       |
| `MorganFingerprint`                  | `'morgan_fingerprint'`              | Morgan fingerprint of radius 3, consisting of 2048 bits packed into a NumPy array. Calculated by RDKit |
| `NumberLipinskiHydrogenDonors`       | `'num_lipinski_hydrogen_donors'`    | Number of Lipinski hydrogen donors, as defined by OpenEye                                              |
| `NumberLipinskiHydrogenAcceptors`    | `'num_lipinski_hydrogen_acceptors'` | Number of Lipinski hydrogen acceptors, as defined by OpenEye                                           |
| `OEXLogP`                            | `'oexlogp'`                         | Modified XLogP calculation as defined by OpenEye                                                       |
| `WildmanCrippenLogP`                 | `'wildmancrippen_logp'`             | Wildman-Crippen LogP calculation as defined by RDKit                                                   |
| `SyntheticAccessibilityScore`        | `'synthetic_accessibility'`         | Synthetic accessibility score, as defined by RDKit                                                     |
| `TwoDimensionalTPSA`                 | `'2d_tpsa'`                         | Two-dimensional Topographical Polar Surface Area, as defined by OpenEye                                |

Molecular properties may be calculated for a SMILES string using:

```python
from chemutils.props import calculate_properties_for_smiles, MolecularWeight, OEXLogP

...

properties = calculate_properties_for_smiles(smiles, {MolecularWeight, OEXLogP})
# properties is a dictionary {'molecular_weight': value, 'oexlogp' :value}

```

They may also be calculated directly from an OpenEye molecule using `calculate_properties_for_oemol`.
