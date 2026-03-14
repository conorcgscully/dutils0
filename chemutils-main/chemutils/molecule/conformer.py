from openeye import oechem, oeomega

OMEGA_RETURN_CODES = {
    oeomega.OEOmegaReturnCode_Success: "Success",
    oeomega.OEOmegaReturnCode_InvalidOptions: "InvalidOptions",
    oeomega.OEOmegaReturnCode_No3DFromCT: "No3DFromCT",
    oeomega.OEOmegaReturnCode_TooManyRotors: "TooManyRotors",
    oeomega.OEOmegaReturnCode_UnspecifiedStereo: "UnspecifiedStereo",
    oeomega.OEOmegaReturnCode_ExplicitHydrogens: "ExplicitHydrogens",
    oeomega.OEOmegaReturnCode_MissingFFParams: "MissingFFParams",
    oeomega.OEOmegaReturnCode_FailedCTBuild: "FailedCTBuild",
    oeomega.OEOmegaReturnCode_FailedTorDrive: "FailedTorDrive",
    oeomega.OEOmegaReturnCode_ExceedsAtomDegree4: "ExceedsAtomDegree4",
    oeomega.OEOmegaReturnCode_FailedDupSetup: "FailedDupSetup",
    oeomega.OEOmegaReturnCode_FailedDGConfGen: "FailedDGConfGen",
    oeomega.OEOmegaReturnCode_FailedTorAssign: "FailedTorAssign",
    oeomega.OEOmegaReturnCode_FailedSmartMatch: "FailedSmartMatch",
    oeomega.OEOmegaReturnCode_FailedFixMatch: "FailedFixMatch",
    oeomega.OEOmegaReturnCode_NoFixedFragment: "NoFixedFragment",
    oeomega.OEOmegaReturnCode_NoValidConfs: "NoValidConfs",
    oeomega.OEOmegaReturnCode_No3DTorDrive: "No3DTorDrive",
    oeomega.OEOmegaReturnCode_Failed: "Failed",
}


class OEOmegaError(RuntimeError):
    pass


def generate_conformers_omega(
    oemol: oechem.OEMolBase, /, *, options: oeomega.OEOmegaOptions | None = None
) -> oechem.OEMCMolBase:
    """
    Run OpenEye's Omega to generate conformers.

    The molecule returned by this function is an `OEMCMol` which contains multiple conformers (coordinates)
    in one OpenEye molecule.

    The returned molecule can have changed, such as gaining hydrogens or atoms being reordered. Passing
    standardised molecules into this will not result in standardised results.

    Args:
        oemol: OpenEye molecule to base conformers on.
        options: Options for OEOmega.

    Returns:
        Multi-conformer OpenEye molecule containing conformers, which can be iterated over using `GetConfs()`.

    Raises:
        OEOmegaError: Error occured when Omega was running.
    """
    options = options or oeomega.OEOmegaOptions()

    omega = oeomega.OEOmega(options)

    oemol_mc = oechem.OEMol(oemol)
    return_code = omega.Build(oemol_mc)

    if return_code == oeomega.OEOmegaReturnCode_Success:
        return oemol_mc
    else:
        reason = OMEGA_RETURN_CODES.get(return_code, f"Unknown Return Code {return_code}")
        raise OEOmegaError(reason)


def generate_single_conformer(
    oemol: oechem.OEMolBase, options: oeomega.OEOmegaOptions | None = None
) -> oechem.OEConfBase:
    """
    Run OpenEye's Omega to generate a single conformer.

    The molecule returned by this function is an `OEConfBase` which represents the single conformer
    (set of coordinates).

    The returned molecule can have changed, such as gaining hydrogens or atoms being reordered. Passing
    standardised molecules into this will not result in standardised results.

    Args:
        oemol: OpenEye molecule to generate a conformer fore.
        options: Options for OEOmega.

    Returns:
        Single-conformer OpenEye molecule.

    Raises:
        OEOmegaError: Error occured when Omega was running.
    """
    options = options or oeomega.OEOmegaOptions()
    options.SetMaxConfs(1)
    (conformer,) = generate_conformers_omega(oemol, options=options).GetConfs()

    return conformer
