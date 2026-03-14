from openeye import oechem


class ReactionError(ValueError):
    pass


class InvalidReactionError(ReactionError):
    pass


class ReactionNotApplicableError(ReactionError):
    pass


class MultipleReactionMatchesError(ReactionError):
    pass


def apply_reaction(
    *,
    reaction: str,
    mol: oechem.OEMolBase,
    correct_valency: bool = False,
    copy_map_idx: bool = False,
) -> oechem.OEMolBase:
    """
    Apply a reaction to a molecule and return the product.

    If the reaction is not applicable to the molecule, or if the reaction can be applied in more
    than one way, an exception is raised.

    Args:
        reaction: Reaction as a SMIRKS string.
        mol: Molecule to apply the reaction to.
        correct_valency: Should the valency of the product be corrected by adding or removing implicit hydrogens.
        copy_map_idx: Should the map indices of the reactants be copied to the product.

    Returns:
        Molecule after the reaction has been applied.

    Raises:
        InvalidReactionError: The reaction is not valid.
        ReactionNotApplicableError: The reaction is not applicable to the molecule.
        MultipleReactionMatchesError: The reaction can be applied in more than one way.
    """
    libgen = oechem.OELibraryGen(reaction)
    if not libgen.IsValid():
        raise InvalidReactionError(f"Invalid reaction `{reaction}`.")

    libgen.SetExplicitHydrogens(False)
    libgen.SetValenceCorrection(correct_valency)
    libgen.SetAssignMapIdx(copy_map_idx)
    libgen.SetStartingMaterial(mol, 0)

    products_iterator = iter(libgen.GetProducts())

    try:
        product = next(products_iterator).CreateCopy()
    except StopIteration:
        raise ReactionNotApplicableError(
            f"Reaction `{reaction}` not applicable to molecule `{oechem.OEMolToSmiles(mol)}`."
        ) from None

    try:
        product = next(products_iterator)
    except StopIteration:
        pass
    else:
        raise MultipleReactionMatchesError(
            f"Reaction `{reaction}` matches more than one part of `{oechem.OEMolToSmiles(mol)}`."
        ) from None

    return product
