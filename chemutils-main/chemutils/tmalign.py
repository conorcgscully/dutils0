import logging
import re
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import tmtools
from openeye import oechem

from chemutils.molecule import transform_molecule
from chemutils.molecule.transform import AffineTransformation
from chemutils.protein.alignment import SequenceAlignment
from chemutils.protein.residues import PROTEIN_CODE_THREE_TO_ONE

STATS_LINE_REGEX = re.compile(
    r"^Aligned length=\s+(?P<aligned_length>\d+), RMSD=\s+(?P<rmsd>\d+.?\d*), Seq_ID=n_identical/n_aligned=\s+(?P<seq_id>\d+.?\d*)$"
)
TMSCORE_LINE_REGEX = re.compile(
    r"^TM-score=\s+(?P<tmscore>\d+.?\d*)\s+\(if normalized by length of Chain_2\)$"
)


logger = logging.getLogger(__name__)


@dataclass
class TMAlignResult:
    aligned_length: int
    rmsd: float
    sequence_identity: float
    tm_score: float
    alignment: SequenceAlignment
    transformation: AffineTransformation


class TMAlignError(RuntimeError):
    pass


def align_with_tmalign(
    *,
    protein: oechem.OEMolBase,
    reference_protein: oechem.OEMolBase,
    chain_id: str | None = None,
    reference_chain_id: str | None = None,
) -> oechem.OEMolBase:
    result = tmalign(
        protein=protein,
        reference_protein=reference_protein,
        chain_id=chain_id,
        reference_chain_id=reference_chain_id,
    )
    return transform_molecule(protein, transformation=result.transformation)


def subset_oemol(
    oemol: oechem.OEMolBase, /, *, subset: oechem.OEUnaryAtomPred, adjust_h_count: bool = False
) -> oechem.OEMolBase:
    new_oemol = oechem.OEGraphMol()
    oechem.OESubsetMol(new_oemol, oemol, subset, adjust_h_count)
    return new_oemol


AMINO_ACIDS = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


def get_protein_and_peptide_chain_ids(oemol: oechem.OEGraphMol, /) -> tuple[set[str], set[str]]:
    """
    Get the list of chain IDs that correspond to protein and peptide chains.

    This function returns two sets of chain IDs: one which contains at least 20
    amino acids (proteins), and one which contains at least 4 amino acids (proteins
    and peptides).
    """
    hv = oechem.OEHierView(oemol)
    amino_acids_per_chain: dict[str, int] = defaultdict(lambda: 0)
    for residue in hv.GetResidues():
        if residue.GetResidueName() in AMINO_ACIDS:
            amino_acids_per_chain[residue.GetOEResidue().GetChainID()] += 1
    return {chain_id for chain_id, count in amino_acids_per_chain.items() if count >= 20}, {
        chain_id for chain_id, count in amino_acids_per_chain.items() if count >= 4
    }


def get_alpha_carbon_coords_and_sequence(
    oemol: oechem.OEMolBase,
) -> tuple[npt.NDArray[np.float32], str]:
    hv = oechem.OEHierView(oemol)
    coords = []
    seq = []
    for residue in hv.GetResidues():
        for atom in residue.GetAtoms():
            if atom.GetName().strip() == "CA" and atom.GetAtomicNum() == 6:
                residue_name = residue.GetResidueName()
                if residue_name not in PROTEIN_CODE_THREE_TO_ONE:
                    continue
                coords.append(oemol.GetCoords(atom))
                seq.append(PROTEIN_CODE_THREE_TO_ONE[residue_name])
    return np.array(coords), "".join(seq)


def tmalign(
    *,
    protein: oechem.OEMolBase,
    reference_protein: oechem.OEMolBase,
    chain_id: str | None = None,
    reference_chain_id: str | None = None,
) -> TMAlignResult:
    """
    Run TMAlign to align two protein chains.

    If either the chain or reference chain have multiple protein chains, then the exact chain must be specified
    using the `chain_id` and `reference_chain_id` arguments. This differs from TMAlign's default behaviour, which
    would just silently take the first chain.

    The transformation matrix in the return object can be used to transform other molecules such as ligands in the
    same manner, using `chemutils.molecule.transform_molecule`.

    Args:
        protein: Protein to align.
        reference_protein: Reference protein to align to.
        chain_id: Chain ID of the protein to align. Not required if there is a single chain in `protein`.
        reference_chain_id: Chain ID of the reference protein to align to. Not required if there is a single chain
            in `reference_protein`.

    Returns:
        Dataclass containing relevant results from `TMAlign`, including:
        * The RMSD of the alignment.
        * The sequence alignment as an `SequenceAlignment` object.
        * The affine transformation matrix required to transform `protein` to be aligned with `reference_protein`.
    """
    # get structures using BioPython

    valid_protein_chain_ids, valid_protein_or_peptide_chain_ids = get_protein_and_peptide_chain_ids(
        protein
    )
    valid_reference_protein_chain_ids, valid_reference_protein_or_peptide_chain_ids = (
        get_protein_and_peptide_chain_ids(reference_protein)
    )

    if chain_id is not None:
        if chain_id not in valid_protein_or_peptide_chain_ids:
            raise TMAlignError(
                f"TMAlign error: Chain ID `{chain_id}` is not a protein or peptide chain."
            )
        protein = subset_oemol(protein, subset=oechem.OEHasChainID(chain_id))

    elif len(valid_protein_chain_ids) > 1:
        raise TMAlignError(
            "TMAlign error: Multiple protein chains found in protein. Specify a chain ID with `chain_id`."
        )

    if reference_chain_id is not None:
        if reference_chain_id not in valid_reference_protein_or_peptide_chain_ids:
            raise TMAlignError(
                f"TMAlign error: Reference chain ID `{reference_chain_id}` is not a protein or peptide chain."
            )
        reference_protein = subset_oemol(
            reference_protein, subset=oechem.OEHasChainID(reference_chain_id)
        )
    elif len(valid_reference_protein_chain_ids) > 1:
        raise TMAlignError(
            "TMAlign error: Multiple protein chains found in reference protein. Specify a chain ID with `reference_chain_id`."
        )

    coords, seq = get_alpha_carbon_coords_and_sequence(protein)
    ref_coords, ref_seq = get_alpha_carbon_coords_and_sequence(reference_protein)

    result = tmtools.tm_align(coords, ref_coords, seq, ref_seq)

    aligned_length = len([ch for ch in result.seqM if ch == ":" or ch == "."])

    num_identical_residues = sum(
        ch1 == ch2 for ch1, ch2 in zip(result.seqxA, result.seqyA, strict=True)
    )

    sequence_identity = num_identical_residues / aligned_length

    transformation = np.eye(4, dtype=np.float32)
    transformation[:3, :3] = result.u
    transformation[:3, 3] = result.t

    return TMAlignResult(
        aligned_length=aligned_length,
        rmsd=result.rmsd,
        sequence_identity=sequence_identity,
        tm_score=result.tm_norm_chain2,
        alignment=SequenceAlignment.from_aligned_sequences(
            aligned_sequence=result.seqxA,
            aligned_reference=result.seqyA,
        ),
        transformation=transformation,
    )
