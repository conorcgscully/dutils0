from .jvm.chemaxon.marvin.calculations import logPPlugin
from .molecule import ConvertibleToCXMol, as_cxmol


def get_logp(mol: ConvertibleToCXMol, /) -> float:
    mol = as_cxmol(mol)

    plugin = logPPlugin()

    plugin.setMolecule(mol)
    plugin.run()

    return plugin.getlogPTrue()  # type: ignore
