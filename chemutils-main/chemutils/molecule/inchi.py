from dataclasses import dataclass

from openeye import oechem


@dataclass
class InChI:
    inchi_key: str
    molformula: str
    connections: str | None = None
    hydrogens: str | None = None
    stereo: str | None = None
    charge: str | None = None
    isotope: str | None = None

    @property
    def inchi(self) -> str:
        return "".join(
            [
                "InChI=1S/",
                self.molformula,
                self.connections if self.connections else "",
                self.hydrogens if self.hydrogens else "",
                self.stereo if self.stereo else "",
                self.charge if self.charge else "",
                self.isotope if self.isotope else "",
            ]
        )


def get_inchi_breakdown(mol: oechem.OEMolBase, /) -> InChI:
    inchi = oechem.OEMolToSTDInChI(mol)
    inchi_key = oechem.OEMolToSTDInChIKey(mol)

    inchi_parts = inchi.split("/")
    assert inchi_parts[0] == "InChI=1S"
    molformula = inchi_parts[1]

    connections = ""
    hydrogens = ""
    stereo = ""
    charge = ""
    isotope = ""

    for part in inchi_parts[2:]:
        if part[0] == "c":
            connections += "/" + part
        elif part[0] == "h":
            hydrogens += "/" + part
        elif part[0] in ["s", "b", "t", "m"]:
            stereo += "/" + part
        elif part[0] in ["q", "p"]:
            charge += "/" + part
        elif part[0] in ["i"]:
            isotope += "/" + part
        else:
            raise ValueError(part)

    assert f"InChI=1S/{molformula}{connections}{hydrogens}{stereo}{charge}{isotope}" == inchi

    result = InChI(
        inchi_key=inchi_key,
        molformula=molformula,
        connections=connections or None,
        hydrogens=hydrogens or None,
        stereo=stereo or None,
        charge=charge or None,
        isotope=isotope or None,
    )
    if result.inchi != inchi:
        raise RuntimeError(f"InChI breakdown was not equal, {result.inchi} != {inchi}")
    return result
