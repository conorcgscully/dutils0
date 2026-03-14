from .jvm.chemaxon.marvin.calculations import logDPlugin
from .molecule import ConvertibleToCXMol, as_cxmol


def get_logd(
    mol: ConvertibleToCXMol, /, *, pH: float = 7.4, consider_tautomers: bool = False
) -> float:
    mol = as_cxmol(mol)

    plugin = logDPlugin()

    plugin.setpH(float(pH))

    if consider_tautomers:
        plugin.setConsiderTautomerization(True)

    plugin.setMolecule(mol)
    plugin.run()

    return plugin.getlogD()  # type: ignore
