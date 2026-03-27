import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md("""
    # Amino acid filtering and clustering notebook

    This notebook loads amino acid building blocks, protects canonical amino acids,
    applies filtering steps, computes fingerprints, and performs greedy centroid
    clustering at a Tanimoto threshold of 0.9.
    """)
    return


@app.cell
def _():
    from ncycle_core import data_utils
    import json
    from pathlib import Path
    import pandas as pd
    import mols2grid

    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator




    from rdkit.Chem import AllChem as Chem
    from rdkit.Chem import PandasTools
    from rdkit.Chem import Descriptors
    from IPython.display import display
    from rdkit.Chem import Draw

    return (
        Chem,
        DataStructs,
        Descriptors,
        Draw,
        MaxMinPicker,
        PandasTools,
        Path,
        data_utils,
        mols2grid,
        pd,
        rdFingerprintGenerator,
    )


@app.cell
def _(Path, data_utils):
    if data_utils.GraphQlDataLoader.is_ready():
        #env variables NCYCLE_STUDIO_API_KEY and NCYCLE_STUDIO_API_URL are properly defined
        data_loader = data_utils.GraphQlDataLoader()
    else:
        #otherwise load them as defined within a json config file
        repo_root = Path(__file__).resolve().parent.parent
        env_config_file = repo_root / "env.prod.json"
        data_loader = data_utils.GraphQlDataLoader(str(env_config_file))

    bbs = data_loader.get_building_blocks_as_json(including_archived=True, including_timestamps=True)
    print(len(bbs))
    return (bbs,)


@app.cell
def _(PandasTools, bbs, pd):
    df_bb = pd.DataFrame.from_records(bbs)
    df_aa = df_bb.query("Type == 'amino_acid'")
    PandasTools.AddMoleculeColumnToFrame(df_aa, smilesCol="Smiles", molCol="ROMol", includeFingerprints=False)
    df_aa.columns
    return df_aa, df_bb


@app.cell
def _(df_bb):
    df_bb.query("Type == 'amino_thiol'")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Arginines
    """)
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Rotatable Bonds
    df_aa['NumRotBonds'] =
    """)
    return


@app.cell
def _(Descriptors, df_aa):
    df_aa["NumRotatableBonds"] = df_aa["ROMol"].apply(
        lambda m: Descriptors.NumRotatableBonds(m) if m is not None else None
    )
    print(df_aa.shape)
    # Filter
    df_aa_1 = df_aa.query("NumRotatableBonds <= 6")
    print(df_aa_1.shape)

    return


@app.cell
def _(Chem, df_aa):

    patt = Chem.MolFromSmarts("[#0]-[N]-[CH3]")

    df_aa["is_N_methyl"] = df_aa["ROMol"].apply(
        lambda s: (s.HasSubstructMatch(patt))
    )
    return


@app.cell
def _(Chem, df_aa):
    _patt = Chem.MolFromSmarts("[#0]-[N]-[CX4;H0]-[C&H3]")

    df_aa["is_a_methyl"] = df_aa["ROMol"].apply(
        lambda s: (s.HasSubstructMatch(_patt))
    )
    # df_aa.query("is_a_methyl == True")
    return


@app.cell
def _(Chem, df_aa):
    _patt = Chem.MolFromSmarts("[#0]-[N]-[CX4]-[C]-[#0]")

    df_aa["is_alpha_aa"] = df_aa["ROMol"].apply(
        lambda s: (s.HasSubstructMatch(_patt))
    )
    df_aa.query("is_alpha_aa == True")
    return


@app.cell
def _(Chem, df_aa):
    _patt = Chem.MolFromSmarts("[#0]-[N]-[CX4]-[C]-[C]-[#0]")

    df_aa["is_beta_aa"] = df_aa["ROMol"].apply(
        lambda s: (s.HasSubstructMatch(_patt))
    )

    return


@app.cell
def _(Chem, df_aa):
    _patt = Chem.MolFromSmarts("[S]-[#0]")

    df_aa["is_S100_aa"] = df_aa["ROMol"].apply(
        lambda s: (s.HasSubstructMatch(_patt))
    )
    return


@app.cell
def _(Chem, Draw, df_aa, mo, pd):
    def _():
        _patt = Chem.MolFromSmarts("a")
        _df = df_aa.query("is_beta_aa == False and is_alpha_aa == False")
        print(len(_df))
        df_arom = df_aa.loc[
            df_aa["ROMol"].apply(lambda mol: mol is not None and mol.HasSubstructMatch(_patt))
        ].copy()
        df_arom_ab = df_arom.query("is_alpha_aa == True or is_beta_aa == True")
        _df = df_arom_ab.copy()
        n = 100
        mols = [Chem.MolFromSmiles(s) for s in _df.head(n)["Smiles"]]
        legends = _df.head(n)["Alias"].fillna("").tolist()

        # Write 4-column CSV matching grid layout
        mols_per_row = 4
        padded = legends + [""] * (-len(legends) % mols_per_row)  # pad to multiple of 4
        rows = [padded[i:i + mols_per_row] for i in range(0, len(padded), mols_per_row)]
        pd.DataFrame(rows, columns=[f"col_{i+1}" for i in range(mols_per_row)]).to_csv(
            "aromatic_aa_grid.csv", index=False
        )
    
        img_svg = Draw.MolsToGridImage(
            mols,
            legends=legends,
            molsPerRow=4,
            subImgSize=(220, 180),
            useSVG=True,
        )
        with open("aromatic_aa_grid.svg", "w") as f:
            f.write(img_svg)

        img_png = Draw.MolsToGridImage(
            mols,
            legends=legends,
            molsPerRow=4,
            subImgSize=(220, 180),
            useSVG=False,
        )
        return mo.md(
            f"""
            Saved to `aromatic_aa_grid.svg`
            {mo.as_html(img_png)}
            """
        )
    _()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Filter
    """)
    return


@app.cell
def _(df_aa):
    df_aa_filtered = (
        df_aa
        .pipe(lambda d: d[d["NumRotatableBonds"] <= 6])
        .pipe(lambda d: d[d["is_N_methyl"].eq(False)])
        .pipe(lambda d: d[d["is_a_methyl"].eq(False)])
        .pipe(lambda d: d[d["is_S100_aa"].eq(False)])
    )
    print(df_aa_filtered.shape)
    return (df_aa_filtered,)


@app.cell
def _(PandasTools, df_aa_filtered):
    PandasTools.AddMurckoToFrame(df_aa_filtered)
    return


@app.cell
def _(PandasTools, df_aa_filtered):
    PandasTools.WriteSDF(df_aa_filtered, "aas_filtered.sdf", molColName='ROMol', properties=df_aa_filtered.columns)
    return


@app.cell
def _(Chem, Draw, df_aa, mo):
    # for IL17
    _patt = Chem.MolFromSmarts("a")

    df_arom = df_aa.loc[
        df_aa["ROMol"].apply(lambda mol: mol is not None and mol.HasSubstructMatch(_patt))
    ].copy()
    df_arom_ab = df_arom.query("is_alpha_aa == True or is_beta_aa == True")


    n=100
    mols = [Chem.MolFromSmiles(s) for s in _df.head(n)["Smiles"]]
    legends = _df.head(n)["Alias"].fillna("").tolist()

    img = Draw.MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=4,
        subImgSize=(220, 180),
        useSVG=False,
    )
    mo.md(
        f"""

        {mo.as_html(img)}
        """
    )
    return (df_arom_ab,)


@app.cell
def _(df_arom_ab):
    for i in df_arom_ab.Alias.values:
        print(f'"{i}",')
    return


@app.cell
def _(Chem):
    mol = """
      -INDIGO-03252608542D

      0  0  0  0  0  0  0  0  0  0  0 V3000
    M  V30 BEGIN CTAB
    M  V30 COUNTS 34 38 1 0 0
    M  V30 BEGIN ATOM
    M  V30 1 O 12.9962 -5.12321 0.0 0
    M  V30 2 C 12.2531 -5.79234 0.0 0
    M  V30 3 O 11.302 -5.48333 0.0 0
    M  V30 4 C 12.461 -6.77049 0.0 0 CFG=2
    M  V30 5 C 11.7179 -7.43962 0.0 0
    M  V30 6 C 10.7668 -7.13061 0.0 0
    M  V30 7 C 10.5589 -6.15246 0.0 0
    M  V30 8 C 9.60783 -5.84344 0.0 0
    M  V30 9 C 8.86469 -6.51257 0.0 0
    M  V30 10 N 7.87017 -6.40804 0.0 0
    M  V30 11 C 7.46343 -7.32159 0.0 0
    M  V30 12 C 8.20658 -7.99072 0.0 0
    M  V30 13 C 9.0726 -7.49072 0.0 0
    M  V30 14 C 10.0237 -7.79974 0.0 0
    M  V30 15 C 8.10205 -8.98524 0.0 0
    M  V30 16 N 7.99752 -9.97976 0.0 0
    M  V30 17 N 13.4121 -7.07951 0.0 0
    M  V30 18 C 14.1552 -6.41038 0.0 0
    M  V30 19 O 13.9473 -5.43223 0.0 0
    M  V30 20 O 15.1063 -6.7194 0.0 0
    M  V30 21 C 15.8494 -6.05027 0.0 0
    M  V30 22 C 16.8005 -6.35928 0.0 0
    M  V30 23 C 17.1095 -7.31034 0.0 0
    M  V30 24 C 16.6095 -8.17636 0.0 0
    M  V30 25 C 17.1095 -9.04239 0.0 0
    M  V30 26 C 18.1095 -9.04239 0.0 0
    M  V30 27 C 18.6095 -8.17636 0.0 0
    M  V30 28 C 18.1095 -7.31034 0.0 0
    M  V30 29 C 18.4185 -6.35928 0.0 0
    M  V30 30 C 19.332 -5.95254 0.0 0
    M  V30 31 C 19.4366 -4.95802 0.0 0
    M  V30 32 C 18.6276 -4.37024 0.0 0
    M  V30 33 C 17.714 -4.77697 0.0 0
    M  V30 34 C 17.6095 -5.7715 0.0 0
    M  V30 END ATOM
    M  V30 BEGIN BOND
    M  V30 1 1 1 2
    M  V30 2 2 2 3
    M  V30 3 1 2 4
    M  V30 4 1 4 5 CFG=1
    M  V30 5 1 5 6
    M  V30 6 4 6 7
    M  V30 7 4 7 8
    M  V30 8 4 8 9
    M  V30 9 4 9 10
    M  V30 10 4 10 11
    M  V30 11 4 11 12
    M  V30 12 4 12 13
    M  V30 13 4 13 9
    M  V30 14 4 13 14
    M  V30 15 4 14 6
    M  V30 16 1 12 15
    M  V30 17 3 15 16
    M  V30 18 1 4 17
    M  V30 19 1 17 18
    M  V30 20 2 18 19
    M  V30 21 1 18 20
    M  V30 22 1 20 21
    M  V30 23 1 21 22
    M  V30 24 1 22 23
    M  V30 25 4 23 24
    M  V30 26 4 24 25
    M  V30 27 4 25 26
    M  V30 28 4 26 27
    M  V30 29 4 27 28
    M  V30 30 4 28 23
    M  V30 31 1 28 29
    M  V30 32 4 29 30
    M  V30 33 4 30 31
    M  V30 34 4 31 32
    M  V30 35 4 32 33
    M  V30 36 4 33 34
    M  V30 37 4 34 29
    M  V30 38 1 34 22
    M  V30 END BOND
    M  V30 BEGIN COLLECTION
    M  V30 MDLV30/STERAC1 ATOMS=(1 4)
    M  V30 END COLLECTION
    M  V30 BEGIN SGROUP
    M  V30 1 DAT 1 ATOMS=(1 10) FIELDNAME=MRV_IMPLICIT_H FIELDDISP="    0.0000   -
    M  V30  0.0000    DA    ALL  1       1  " FIELDDATA="IMPL_H1"
    M  V30 END SGROUP
    M  V30 END CTAB
    M  END
    """
    Chem.MolToCXSmiles(Chem.MolFromMolBlock(mol))
    return


@app.cell(hide_code=True)
def _(DataStructs, MaxMinPicker, df_aa_filtered, pd, rdFingerprintGenerator):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    def _():

        def mol_to_fp(mol):
            return mfpgen.GetFingerprint(mol)
    
    
        def diverse_order_for_group(group: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
            fps = group["fp"].tolist()
            n = len(fps)
    
            if n <= 1:
                return group.copy()
    
            picker = MaxMinPicker()
    
            def dist(i: int, j: int) -> float:
                return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    
            order = list(picker.LazyPick(dist, n, n, seed=seed))
            return group.iloc[order].copy()
    
    
        df = df_aa_filtered.copy()
        df = df[df["ROMol"].notna() & df["Murcko_SMILES"].notna()].copy()
        df["fp"] = df["ROMol"].map(mol_to_fp)
    
        picker = MaxMinPicker()
    
        # 1 pick per scaffold first
        first_picks = []
        remaining = []
    
        for scaffold, g in df.groupby("Murcko_SMILES", sort=False):
            ordered = diverse_order_for_group(g, seed=0)
            first_picks.append(ordered.iloc[0])
            if len(ordered) > 1:
                remaining.append(ordered.iloc[1:])
    
        df_first = pd.DataFrame(first_picks)
        df_remaining = (
            pd.concat(remaining, axis=0).reset_index(drop=True)
            if remaining
            else pd.DataFrame(columns=df.columns)
        )
    
        # case 1: enough scaffolds already
        if len(df_first) >= 300:
            fps = df_first["fp"].tolist()
            n_fps = len(fps)
    
            def distij(i: int, j: int) -> float:
                return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    
            picks = list(picker.LazyPick(distij, n_fps, 300, seed=0))
            df_selected = df_first.iloc[picks].reset_index(drop=True)
    
        # case 2: need more molecules
        else:
            pool = pd.concat([df_first, df_remaining], axis=0).reset_index(drop=True)
            fps = pool["fp"].tolist()
            n_pool = len(pool)
    
            # simple fallback: take scaffold reps first, then add diverse picks from remainder
            already_n = len(df_first)
            remaining_n = min(300 - already_n, len(df_remaining))
    
            if remaining_n > 0:
                rem_fps = df_remaining["fp"].tolist()
    
                def dist_rem(i: int, j: int) -> float:
                    return 1.0 - DataStructs.TanimotoSimilarity(rem_fps[i], rem_fps[j])
    
                rem_picks = list(picker.LazyPick(dist_rem, len(rem_fps), remaining_n, seed=0))
                df_extra = df_remaining.iloc[rem_picks]
                df_selected = pd.concat([df_first, df_extra], axis=0).reset_index(drop=True)
            else:
                df_selected = df_first.reset_index(drop=True)
    
        return df_selected
    _()
    return


@app.cell
def _(df_aa_filtered):
    df_aa_filtered.columns
    return


@app.cell
def _(df_aa, mols2grid):
    mols2grid.display(
        df_aa.head(),
        mol_col="ROMol",
        subset=["img", "Alias"],
    )
    return


@app.cell
def _(df_aa):
    df_aa.query("is_beta_aa == False and is_alpha_aa == False ")
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Protect canonical AAs
    """)
    return


@app.cell
def _(df_aa):
    canon_list = [
        "Ala",
        "Asp",
        "Glu",
        "Phe",
        "Gly",
        "His",
        "Ile",
        "Lys",
        "Leu",
        "Asn",
        "Pro",
        "Gln",
        "Arg(Pbf)",
        "Ser",
        "Thr",
        "Val",
        "Trp",
        "Tyr",
        "D-Ala",
        "D-Asp",
        "D-Glu",
        "D-Phe",
        "D-His",
        "D-Ile",
        "D-Lys",
        "D-Leu",
        "D-Asn",
        "D-Pro",
        "D-Gln",
        "D-Arg",
        "D-Ser",
        "D-Thr",
        "D-Val",
        "D-Trp",
        "D-Tyr"
    ]
    for _a in canon_list:
        if _a not in df_aa['Alias'].values:
            print(_a)
    return


@app.cell
def _(PandasTools, df_raw):
    smiles_col = _resolve_smiles_column(df_raw)
    df_with_mols = df_raw.copy()
    PandasTools.AddMoleculeColumnToFrame(df_with_mols, smilesCol=smiles_col, molCol="ROMol", includeFingerprints=False)
    dropped_mol_count = int(df_with_mols["ROMol"].isna().sum())
    df_valid = df_with_mols.loc[df_with_mols["ROMol"].notna()].copy()
    assert "ROMol" in df_valid.columns
    assert len(df_valid) > 0
    return (df_valid,)


@app.cell
def _(mo):
    mo.md("""
    ## 3) Protect canonical amino acids

    Build a frozen canonical-SMILES priority set for the 20 L-amino acids plus
    their D-enantiomers; these rows bypass the early filters.
    """)
    return


@app.cell
def _(Chem, df_valid):
    protected = _protected_smiles_set()
    df_tagged = df_valid.copy()
    df_tagged["canonical_smiles"] = df_tagged["ROMol"].map(lambda m: Chem.MolToSmiles(m, isomericSmiles=True, canonical=True))
    df_tagged["is_protected"] = df_tagged["canonical_smiles"].isin(protected)
    assert "is_protected" in df_tagged.columns
    return (df_tagged,)


@app.cell
def _(mo):
    mo.md("""
    ## 4) Category filter

    Retain non-protected rows where `category` matches `amino_acid`
    after lowercasing and trimming whitespace.
    """)
    return


@app.cell
def _(df_tagged):
    df_protected = df_tagged.loc[df_tagged["is_protected"]].copy()
    df_non_protected = df_tagged.loc[~df_tagged["is_protected"]].copy()
    df_after_category = _filter_category(df_non_protected)
    assert len(df_after_category) >= 0
    return df_after_category, df_protected


@app.cell
def _(df_after_category, mo):
    mo.md(f"""
    Rows after category filter (non-protected only): **{len(df_after_category)}**
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    ## 5) SMARTS check for N-methyl amino acids

    The SMARTS flag is computed and retained for interpretation. The inclusion logic
    intentionally keeps both N-methyl and non-N-methyl amino acids.
    """)
    return


@app.cell
def _(config, df_after_category):
    df_after_smarts = _annotate_n_methyl(df_after_category, config.n_methyl_smarts)
    assert "is_n_methyl" in df_after_smarts.columns
    return (df_after_smarts,)


@app.cell
def _(df_after_smarts, mo):
    n_methyl_count = int(df_after_smarts["is_n_methyl"].sum()) if len(df_after_smarts) else 0
    mo.md(
        f"""
        N-methyl matches in category-filtered rows: **{n_methyl_count}**.

        SMARTS filter keeps both matching and non-matching rows by design.
        """
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## 6) Rotatable bond filter

    Keep non-protected rows with rotatable bonds <= configured limit.
    """)
    return


@app.cell
def _(config, df_after_smarts):
    before_rotb = len(df_after_smarts)
    df_after_rotb = _filter_rotatable_bonds(df_after_smarts, config.max_rot_bonds)
    after_rotb = len(df_after_rotb)
    assert len(df_after_rotb) >= 0
    return after_rotb, before_rotb, df_after_rotb


@app.cell
def _(after_rotb, before_rotb, mo):
    mo.md(f"""
    Rotatable bond filter: **{before_rotb} -> {after_rotb}**
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    ## 7) Merge protected and filtered sets, then deduplicate

    Combine protected rows with filtered non-protected rows and deduplicate by
    canonical SMILES.
    """)
    return


@app.cell
def _(df_after_rotb, df_protected, pd):
    df_merged = pd.concat([df_protected, df_after_rotb], axis=0, ignore_index=True).copy()
    before_dedup = len(df_merged)
    df_dedup = df_merged.drop_duplicates(subset=["canonical_smiles"], keep="first").copy()
    after_dedup = len(df_dedup)
    assert len(df_dedup) > 0
    return after_dedup, before_dedup, df_dedup


@app.cell
def _(after_dedup, before_dedup, mo):
    mo.md(f"""
    Deduplication by canonical SMILES: **{before_dedup} -> {after_dedup}**
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    ## 8) Fingerprints

    Compute Morgan fingerprints for each molecule.
    """)
    return


@app.cell
def _(config, df_dedup):
    df_with_fp = _compute_fingerprints(df_dedup, config.fp_radius, config.fp_nbits)
    assert "fp" in df_with_fp.columns
    return (df_with_fp,)


@app.cell
def _(mo):
    mo.md("""
    ## 9) Greedy centroid clustering (Tanimoto)

    Protected rows are ordered first; each molecule joins an existing centroid
    cluster if similarity >= threshold, otherwise it starts a new cluster.
    """)
    return


@app.cell
def _(config, df_with_fp, greedy_centroid_cluster):
    centroid_df, cluster_sizes = greedy_centroid_cluster(df_with_fp, config.tanimoto_threshold)
    assert len(centroid_df) > 0
    assert "cluster_id" in centroid_df.columns
    assert "cluster_size" in centroid_df.columns
    return (centroid_df,)


@app.cell
def _(mo):
    mo.md("""
    ## 10) Output and visualization

    Show summary counts and attempt a grid image render for centroid molecules.
    """)
    return


@app.cell
def _(
    after_dedup,
    after_rotb,
    centroid_df,
    df_after_category,
    df_after_smarts,
    df_raw,
    np,
    pd,
):
    summary_df = pd.DataFrame(
        [
            {"step": "total_input", "count": int(len(df_raw))},
            {"step": "after_category_filter", "count": int(len(df_after_category))},
            {"step": "after_smarts_filter", "count": int(len(df_after_smarts))},
            {"step": "after_rotb_filter", "count": int(after_rotb)},
            {"step": "after_dedup", "count": int(after_dedup)},
            {"step": "final_centroids", "count": int(len(centroid_df))},
        ]
    )
    # Uses numpy for a simple aggregate metric exposed in notebook output.
    summary_total = int(np.sum(summary_df["count"].values))
    return summary_df, summary_total


@app.cell
def _(centroid_df, mo, summary_df, summary_total):
    grid_html = _grid_image_html(centroid_df, mol_col="ROMol", legends_col="cluster_id")
    summary_view = mo.plain(summary_df)
    if grid_html is not None:
        viz_view = mo.Html(grid_html)
    else:
        viz_view = mo.plain(centroid_df)
    mo.vstack(
        [
            mo.md(f"Summary count checksum (sanity metric): **{summary_total}**"),
            summary_view,
            viz_view,
            mo.plain(centroid_df.head(20)),
        ]
    )
    return


@app.cell
def _(centroid_df, input_path, mo):
    mo.md(
        """
        ### Optional export

        Run the next cell manually to write centroid records to JSON.
        """
    )
    export_default_path = input_path.with_name(f"{input_path.stem}.centroids.json")
    mo.plain({"export_path": str(export_default_path), "rows": len(centroid_df)})
    return


@app.cell
def _(mo):
    # FIXME: This cell is intentionally side-effectful; keep commented unless export is desired.
    # export_default_path.write_text(
    #     json.dumps(_records_for_export(centroid_df), indent=2),
    #     encoding="utf-8",
    # )
    mo.md("""
    Export is disabled by default. Uncomment in-cell code to write JSON.
    """)
    return


if __name__ == "__main__":
    app.run()
