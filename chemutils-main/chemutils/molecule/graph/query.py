"""Convert `AtomProperty` and `BondProperty` objects into query expressions that can be used by OpenEye."""

from openeye import oechem

from ..equal import AtomProperty, BondProperty


def get_query_from_oemol(
    query: oechem.OEMolBase,
    /,
    *,
    atom_properties: AtomProperty,
    bond_properties: BondProperty,
) -> oechem.OEQMol:
    """Create an OpenEye query molecule."""
    atomexpr = get_query_for_atom_properties(atom_properties)
    bondexpr = get_query_for_bond_properties(bond_properties)

    query_mol = oechem.OEQMol(query)
    query_mol.BuildExpressions(atomexpr, bondexpr)

    return query_mol


def get_query_for_atom_properties(properties: AtomProperty, /) -> int:
    """Convert an `AtomProperty` object into a query expression for atoms."""
    atomexpr = 0
    for property in [prop for prop in AtomProperty if prop & properties]:
        atomexpr |= _get_query_for_atom_property(property)
    return atomexpr


def get_query_for_bond_properties(properties: BondProperty, /) -> int:
    """Convert a `BondProperty` object into a query expression for bonds."""
    bondexpr = 0
    for property in [prop for prop in BondProperty if prop & properties]:
        bondexpr |= _get_query_for_bond_property(property)
    return bondexpr


def _get_query_for_atom_property(property: AtomProperty, /) -> int:
    match property:
        case AtomProperty.Aromatic:
            return oechem.OEExprOpts_Aromaticity  # type: ignore
        case AtomProperty.AtomicNumber:
            return oechem.OEExprOpts_AtomicNumber  # type: ignore
        case AtomProperty.Stereo:
            return oechem.OEExprOpts_Chiral  # type: ignore
        case AtomProperty.Degree:
            return oechem.OEExprOpts_Degree  # type: ignore
        case AtomProperty.FormalCharge:
            return oechem.OEExprOpts_StrictFormalCharge  # type: ignore
        case AtomProperty.HCount:
            return oechem.OEExprOpts_HCount  # type: ignore
        case AtomProperty.HeavyDegree:
            return oechem.OEExprOpts_HvyDegree  # type: ignore
        case AtomProperty.Hybridisation:
            return oechem.OEExprOpts_Hybridization  # type: ignore
        case AtomProperty.ImplicitHCount:
            return oechem.OEExprOpts_ImplicitHCount  # type: ignore
        case AtomProperty.Isotope:
            return oechem.OEExprOpts_Mass  # type: ignore
        case AtomProperty.IsInRing:
            return oechem.OEExprOpts_RingMember  # type: ignore
        case AtomProperty.Valence:
            return oechem.OEExprOpts_Valence  # type: ignore
        case _:
            raise ValueError(f"Cannot express property `{property}` as an `atomexpr` query.")


def _get_query_for_bond_property(property: BondProperty, /) -> int:
    match property:
        case BondProperty.Aromatic:
            return oechem.OEExprOpts_Aromaticity  # type: ignore
        case BondProperty.Order:
            return oechem.OEExprOpts_BondOrder  # type: ignore
        case BondProperty.Stereo:
            return oechem.OEExprOpts_Chiral  # type: ignore
        case BondProperty.IsInRing:
            return oechem.OEExprOpts_RingMember  # type: ignore
        case _:
            raise ValueError(f"Cannot express property `{property}` as a `bondexpr` query.")
