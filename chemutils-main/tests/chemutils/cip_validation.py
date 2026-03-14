from pathlib import Path

from openeye import oechem

from chemutils import fs

CIP_VALIDATION_OEMOLS = []
CIP_VALIDATION_TETRA_STEREO_CASES = []

CIP_VALIDATION_OEMOLS_3D = []

for oemol in fs.read_molecules(Path(__file__).parent / "cip_validation_2d.sdf"):
    CIP_VALIDATION_OEMOLS.append(oemol)
    rules = oechem.OEGetSDData(oemol, "SEQUENCE_RULES")
    if "4a" in rules or "4b" in rules or "4c" in rules or "5" in rules or "6" in rules:
        continue
    for label in oechem.OEGetSDData(oemol, "CIP_LABELS").split():
        cip_label = label[-1]
        if cip_label in ("R", "S"):
            atom_idx = int(label[:-1]) - 1
            if oemol.GetTitle() in ("VS031", "VS032"):
                cip_label = "R"
            CIP_VALIDATION_TETRA_STEREO_CASES.append((oemol, atom_idx, cip_label))

for oemol in fs.read_molecules(Path(__file__).parent / "cip_validation_3d.sdf"):
    CIP_VALIDATION_OEMOLS_3D.append(oemol)
