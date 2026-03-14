from collections.abc import Callable, Iterable
from typing import TypeAlias

from openeye import oechem

from .evaluate import evaluate_selection

AtomsSelection: TypeAlias = (
    int
    | oechem.OEAtomBase
    | oechem.OEUnaryAtomBoolFunc
    | oechem.OEAtomBondSet
    | oechem.OEAtomIter
    | oechem.OEMatchBase
    | Callable[[oechem.OEAtomBase], bool]
    | Iterable[int | oechem.OEAtomBase]
    | str
)


def _selection_as_atombondset_or_unaryatompred(
    selection: AtomsSelection, /, *, mol: oechem.OEMolBase
) -> oechem.OEAtomBondSet | oechem.OEUnaryAtomPred:
    match selection:
        case int():
            return oechem.OEHasAtomIdx(selection)
        case oechem.OEAtomBase():
            new_selection = oechem.OEAtomBondSet()
            new_selection.AddAtom(selection)
            return new_selection
        case str():
            return evaluate_selection(oemol=mol, query=selection)
        case oechem.OEUnaryAtomBoolFunc():
            return selection
        case oechem.OEAtomBondSet():
            return selection
        case iterable if isinstance(iterable, Iterable):
            new_selection = oechem.OEAtomBondSet()
            for item in iterable:
                if isinstance(item, int):
                    new_selection.AddAtom(mol.GetAtom(oechem.OEHasAtomIdx(item)))
                else:
                    new_selection.AddAtom(item)
            return new_selection
        case _ if callable(selection):
            return oechem.PyAtomPredicate(selection)
        case _:
            raise TypeError(f"Invalid selection {selection}")


def _selection_as_atombondset(
    selection: AtomsSelection, /, *, mol: oechem.OEMolBase
) -> oechem.OEAtomBondSet:
    atombondset_or_unaryatompred = _selection_as_atombondset_or_unaryatompred(selection, mol=mol)
    if isinstance(atombondset_or_unaryatompred, oechem.OEAtomBondSet):
        return atombondset_or_unaryatompred
    return _selection_as_atombondset(
        [atom for atom in mol.GetAtoms() if atombondset_or_unaryatompred(atom)],
        mol=mol,
    )
