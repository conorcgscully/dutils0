from ..struc import Molecule

class MolImporter:
    @classmethod
    def importMol(self, s: str, opts: str, /) -> Molecule: ...

class MolExporter:
    @classmethod
    def exportToFormat(self, molecule: Molecule, fmt: str, /) -> str: ...
