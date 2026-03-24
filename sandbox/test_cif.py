import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")


@app.cell
def _():
    from dutils0.mol.read import read_protein_bytes
    from dutils0.mol.write import write_protein_str
    from dutils0.mol.convert import protein_to_rdkit, gemmi_from_rdkit
    from pymol import cmd
    import os
    import gemmi
    from rdkit import Chem
    import tempfile
    import os
    from rdsl import select_atoms, select_atom_ids, select_molecule, SelectionResult
    from ipymolstar import PDBeMolstar

    from rdkit.Chem import Draw
    import marimo as mo

    def rodraw(mol, size=(400, 300)):
        return mo.image(Draw.MolToImage(mol, size=size))

    st = read_protein_bytes('dutils0/tests/data/il5r_boltz_model.cif')
    st = read_protein_bytes('dutils0/tests/data/6H41.cif')
    # pymol_obj = gemmi_to_pymol(st, "m1")
    return Chem, PDBeMolstar, protein_to_rdkit, rodraw, select_molecule, st


@app.cell
def _(protein_to_rdkit, st):
    # _o1 = pymol_select(pymol_obj, "resi 200", "res200")
    pr_romol = protein_to_rdkit(st)
    pr_romol.GetNumAtoms()
    return (pr_romol,)


@app.cell
def _(Chem, PDBeMolstar, pr_romol):
    molblock = Chem.MolToPDBBlock(pr_romol)  # mol is an RDKit Mol with a conformer

    view = PDBeMolstar(
        custom_data={
            "data": molblock,
            "format": "pdb",
            "binary": False,
        },
        theme="light",
    )

    view
    return


@app.cell
def _(pr_romol, rodraw, select_molecule):

    # Chem.SanitizeMol(pr_romol)
    result = select_molecule(pr_romol, "chain A and resi 210")
    rodraw(result.mol)
    return


@app.cell
def _(Chem, rodraw):
    _m = Chem.MolFromSmiles('c1cnccc1O')
    for _a in _m.GetAtoms():
        print(dir(_a))
    rodraw(_m)
    return


@app.cell
def _(Chem):

    _m = Chem.MolFromSmiles('c1cnccc1O')

    mol = Chem.AddHs(_m)
    # Set atom names
    canonical_order = Chem.CanonicalRankAtoms(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    for atom, can_idx in zip(mol.GetAtoms(), canonical_order):
        atom_name = atom.GetSymbol().upper() + str(can_idx + 1)
        atom.SetProp("name", atom_name)
    return (canonical_order,)


@app.cell
def _(canonical_order):
    list(canonical_order)
    return


if __name__ == "__main__":
    app.run()
