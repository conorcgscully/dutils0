import fsutils as fs
from tqdm import tqdm

from chemutils.props import (
    ALL_PROPERTIES,
    calculate_properties_for_representative_smiles,
)

final = []

for drug in tqdm(fs.read_yaml("tests/data/drugs.yaml")):
    result = calculate_properties_for_representative_smiles(
        drug["smiles"], properties=ALL_PROPERTIES
    )
    result["morgan_fingerprint"] = "".join(f"{x:02x}" for x in result["morgan_fingerprint"])
    result = {"name": drug["name"]} | result
    for key in result:
        if isinstance(result[key], float):
            result[key] = round(result[key], 2)
    final.append(result)

fs.write_yaml("tests/data/drugs.yaml", contents=final, sort_keys=True)
