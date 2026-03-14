from itertools import groupby
from typing import Any, cast, overload

from openeye import oechem


@overload
def get_sd_tags(oemol: oechem.OEGraphMol, /, name: str) -> list[str]: ...


@overload
def get_sd_tags(oemol: oechem.OEGraphMol, /, name: None) -> dict[str, str | list[str]]: ...


def get_sd_tags(
    oemol: oechem.OEGraphMol, /, name: str | None = None
) -> dict[str, str | list[str]] | list[str]:
    """
    Get all the SD tags for a molecule.

    If `name` is provided, returns a list of values for the tag with the given name.

    If `name` is not provided, returns a dictionary of all tags with their values. If more than one tag with the same
    name exists, the result is a list of its values.

    Raises:
        LookupError: If no tags with the given name are found.
    """
    results = [
        (data_pair.GetTag(), data_pair.GetValue()) for data_pair in oechem.OEGetSDDataPairs(oemol)
    ]
    if name is not None:
        lst = [value for tag, value in results if tag == name]
        if len(lst) == 0:
            raise LookupError(f"Getting SD tags failed: no tags with name `{name}`.")
        return lst
    else:
        results = sorted(results)
        dct = {}
        for tag, values_group in groupby(results, lambda pair: pair[0]):
            values = list(values_group)
            if len(values) == 1:
                dct[tag] = values[0][1]
            else:
                dct[tag] = [value for _, value in values]
        return dct


def has_sd_tag(oemol: oechem.OEGraphMol, /, name: str) -> bool:
    """Does the given molecule have at least one SD tag with the given name."""
    return cast(bool, oechem.OEHasSDData(oemol, name))


def get_sd_tag(oemol: oechem.OEGraphMol, /, name: str) -> str:
    """
    Get the SD tag with the given name for a molecule.

    Raises:
        LookupError: If no tags with the given name are found.
        LookupError: If multiple tags with the given name are found.
    """
    if not has_sd_tag(oemol, name=name):
        raise LookupError(f"Getting SD tag failed: no tags with name `{name}`.")
    values = get_sd_tags(oemol, name=name)
    if len(values) > 1:
        raise LookupError(f"Getting SD tag failed: multiple tags with name `{name}`.")
    return values[0]


def set_sd_tag(oemol: oechem.OEGraphMol, /, *, name: str, value: Any) -> None:
    if has_sd_tag(oemol, name=name) and len(get_sd_tags(oemol, name=name)) > 1:
        raise ValueError(
            f"Setting SD tag failed: multiple tags with name `{name}`. Use `append_sd_tag` to append a new value."
        )
    oechem.OESetSDData(oemol, name, str(value))


def append_sd_tag(oemol: oechem.OEGraphMol, /, *, name: str, value: str) -> None:
    oechem.OEAddSDData(oemol, name, value)


def delete_sd_tags(oemol: oechem.OEGraphMol, /, name: str | None = None) -> None:
    if name is not None:
        oechem.OEDeleteSDData(oemol, name)
    else:
        oechem.OEClearSDData(oemol)
