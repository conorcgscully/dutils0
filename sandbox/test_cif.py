import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")


@app.cell
def _():
    from dutils0.rough.read import read_protein_bytes
    from dutils0.rough.write import write_protein_str
    from dutils0.rough.convert import protein_to_rdkit, gemmi_from_rdkit
    from pymol import cmd
    import os
    import gemmi
    from rdkit import Chem
    import tempfile
    import os
    st = read_protein_bytes('dutils0/tests/data/6H41.cif')
    # pymol_obj = gemmi_to_pymol(st, "m1")
    return gemmi_from_rdkit, protein_to_rdkit, st, write_protein_str


@app.cell
def _(protein_to_rdkit, st):
    # _o1 = pymol_select(pymol_obj, "resi 200", "res200")
    pr_romol = protein_to_rdkit(st)
    pr_romol.GetNumAtoms()
    return (pr_romol,)


@app.cell
def _(gemmi_from_rdkit, pr_romol, write_protein_str):
    pr = gemmi_from_rdkit(pr_romol)
    write_protein_str(pr, fmt="cif")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
