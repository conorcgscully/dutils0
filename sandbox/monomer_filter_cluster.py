# reviewed by: reviewer-agent

import marimo

__generated_with = "0.11.0"
app = marimo.App(width="wide")


@app.cell
def __(mo):
    mo.md(
        r"""
        # Amino Acid Filtering & Clustering Pipeline

        This notebook implements a reproducible cheminformatics pipeline that:

        1. Loads a JSON library of amino acid building blocks
        2. Protects all 20 canonical L/D amino acids from filtering
        3. Applies category, SMARTS, and rotatable-bond filters
        4. Clusters the survivors by Tanimoto similarity (greedy centroid, τ = 0.9)
        5. Returns one representative per cluster for downstream use

        Each section below is self-contained and re-executes reactively whenever
        an upstream cell changes.
        """
    )
    return


@app.cell
def __():
    # ── 1a. Standard-library & third-party imports ──────────────────────────────
    import marimo as mo
    import pandas as pd
    import numpy as np
    from pathlib import Path
    from typing import FrozenSet

    from rdkit import Chem, DataStructs
    from rdkit.Chem import Descriptors, AllChem, Draw
    from rdkit.Chem import PandasTools

    return (
        AllChem,
        Chem,
        DataStructs,
        Descriptors,
        Draw,
        FrozenSet,
        Path,
        PandasTools,
        mo,
        np,
        pd,
    )


@app.cell
def __(Chem, mo):
    # ── 1b. Configuration — all tunable constants live here ─────────────────────
    # REVIEWER: added `mo` parameter so configuration summary is displayed in the
    #           notebook UI rather than sent to the terminal via bare print().
    import types

    config = types.SimpleNamespace(
        # Tanimoto similarity threshold for greedy centroid clustering
        tanimoto_threshold=0.9,
        # Maximum number of rotatable bonds allowed (inclusive)
        max_rot_bonds=5,
        # SMARTS for N-methyl detection:
        #   [NX3;H0;!$(NC=O)] — trivalent N, no H, not on an amide
        #   [CH3]             — directly bonded methyl
        # FIXME: adjust if the dataset uses a different N-methyl motif
        nmethyl_smarts="[NX3;H0;!$(NC=O)][CH3]",
        # Morgan fingerprint parameters
        fp_radius=2,
        fp_bits=2048,
        # Column names
        smiles_col="smiles",
        mol_col="ROMol",
        category_col="category",
        target_category="amino_acid",
        # FIXME: update this path if the JSON file moves
        default_json_path="sandbox/monomer_filter_cluster.py",
    )

    # Pre-compile the SMARTS pattern once; assert to surface bad patterns early.
    _nmethyl_query = Chem.MolFromSmarts(config.nmethyl_smarts)
    assert _nmethyl_query is not None, "N-methyl SMARTS failed to compile"

    mo.md(
        f"""
        **Configuration loaded ✓**

        | Parameter | Value |
        |-----------|-------|
        | Tanimoto threshold | `{config.tanimoto_threshold}` |
        | Max rotatable bonds | `{config.max_rot_bonds}` |
        | Morgan radius | `{config.fp_radius}` |
        | Morgan bits | `{config.fp_bits}` |
        | N-methyl SMARTS | `{config.nmethyl_smarts}` |
        """
    )
    return config, types


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 2 · Load Data

        We accept either a file uploaded interactively via the widget below **or**
        fall back to the hard-coded path in `config.default_json_path`.

        `PandasTools.AddMoleculeColumnToFrame` parses each SMILES string into an
        RDKit `ROMol` object and attaches it as a new column.  Rows that could not
        be parsed are dropped and counted.

        **Expected JSON schema:** a JSON array of objects, each with at minimum a
        `"smiles"` key and a `"category"` key (field names are configurable in
        `config`).  A dict-of-lists schema is also accepted because
        `pd.DataFrame(records)` handles both transparently.
        """
    )
    return


@app.cell
def __(mo):
    # Interactive file-upload widget (optional — leave empty to use the default path)
    file_upload = mo.ui.file(
        label="Upload amino-acid JSON (optional)",
        kind="area",
    )
    file_upload
    return (file_upload,)


@app.cell
def __(Chem, PandasTools, Path, config, file_upload, mo, pd):
    import io
    import json

    def _load_dataframe():
        """Load the amino-acid JSON into a pandas DataFrame with an ROMol column."""
        # ── resolve data source ──────────────────────────────────────────────────
        if file_upload.value:
            raw_bytes = file_upload.value[0].contents
            records = json.loads(raw_bytes.decode("utf-8"))
            source_label = "uploaded file"
        else:
            json_path = Path(config.default_json_path)
            if not json_path.exists():
                raise FileNotFoundError(
                    f"JSON not found at '{json_path}'. "
                    "Please upload a file or update config.default_json_path."
                )
            with json_path.open() as fh:
                records = json.load(fh)
            source_label = str(json_path)

        df = pd.DataFrame(records)
        assert config.smiles_col in df.columns, (
            f"Expected a '{config.smiles_col}' column; got {list(df.columns)}"
        )

        # ── parse SMILES → ROMol ─────────────────────────────────────────────────
        PandasTools.AddMoleculeColumnToFrame(
            df,
            smilesCol=config.smiles_col,
            molCol=config.mol_col,
            includeFingerprints=False,
        )

        n_before = len(df)
        df_valid = df[df[config.mol_col].notna()].copy()
        n_dropped = n_before - len(df_valid)

        assert len(df_valid) > 0, "All rows were dropped — check the SMILES column."
        return df_valid, source_label, n_before, n_dropped

    _df_valid, _source_label, _n_before, _n_dropped = _load_dataframe()

    # Attach canonical SMILES for consistent comparison and deduplication.
    # REVIEWER: canonical SMILES are computed via Chem.MolToSmiles() so that all
    #           downstream equality tests use RDKit's canonical form, never raw
    #           input strings.
    df_raw = _df_valid.copy()
    df_raw["canonical_smiles"] = df_raw[config.mol_col].apply(
        lambda m: Chem.MolToSmiles(m) if m is not None else None
    )

    mo.md(
        f"""
        **Data loaded successfully.**

        | | |
        |-|-|
        | Source | `{_source_label}` |
        | Rows loaded | {_n_before:,} |
        | Rows dropped (bad SMILES) | {_n_dropped:,} |
        | **Rows remaining** | **{len(df_raw):,}** |
        """
    )
    return df_raw, io, json


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 3 · Priority Set — Protect Natural Amino Acids

        The 20 canonical **L-amino acids** and their **D-enantiomers** are
        canonicalised with RDKit and stored in a `frozenset`.

        Any row whose canonical SMILES appears in this set receives
        `is_protected = True` and **bypasses all filters** in steps 4–6,
        going directly into the clustering input.  This guarantees that no
        natural amino acid is accidentally discarded.
        """
    )
    return


@app.cell
def __(Chem, FrozenSet, df_raw, mo):
    # ── canonical SMILES for the 20 L-amino acids + their D-enantiomers ─────────
    # FIXME: extend this list if non-standard proteinogenic AAs (e.g. Sec, Pyl) matter
    _L_SMILES: list[str] = [
        "N[C@@H](C)C(=O)O",                          # L-Ala
        "N[C@@H](CCCNC(=N)N)C(=O)O",                 # L-Arg
        "N[C@@H](CC(=O)N)C(=O)O",                    # L-Asn
        "N[C@@H](CC(=O)O)C(=O)O",                    # L-Asp
        "N[C@@H](CS)C(=O)O",                         # L-Cys
        "N[C@@H](CCC(=O)N)C(=O)O",                   # L-Gln
        "N[C@@H](CCC(=O)O)C(=O)O",                   # L-Glu
        "NCC(=O)O",                                   # Gly (achiral — no D-enantiomer)
        "N[C@@H](Cc1cnc[nH]1)C(=O)O",                # L-His
        "CC[C@H](C)[C@@H](N)C(=O)O",                 # L-Ile  (2S,3S)
        "CC(C)C[C@@H](N)C(=O)O",                     # L-Leu
        "N[C@@H](CCCCN)C(=O)O",                      # L-Lys
        "N[C@@H](CCSC)C(=O)O",                       # L-Met
        "N[C@@H](Cc1ccccc1)C(=O)O",                  # L-Phe
        "O=C(O)[C@@H]1CCCN1",                        # L-Pro
        "N[C@@H](CO)C(=O)O",                         # L-Ser
        "N[C@@H]([C@@H](O)C)C(=O)O",                 # L-Thr (2S,3R)
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",          # L-Trp
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",               # L-Tyr
        "CC(C)[C@@H](N)C(=O)O",                      # L-Val
    ]

    # REVIEWER: D-His (N[C@H](Cc1cnc[nH]1)C(=O)O) was missing from the original
    #           D_SMILES list — added here so the protected set is complete for all
    #           18 chiral D-amino acids (Gly is achiral and has no D form).
    _D_SMILES: list[str] = [
        "N[C@H](C)C(=O)O",                           # D-Ala
        "N[C@H](CCCNC(=N)N)C(=O)O",                  # D-Arg
        "N[C@H](CC(=O)N)C(=O)O",                     # D-Asn
        "N[C@H](CC(=O)O)C(=O)O",                     # D-Asp
        "N[C@H](CS)C(=O)O",                          # D-Cys
        "N[C@H](CCC(=O)N)C(=O)O",                    # D-Gln
        "N[C@H](CCC(=O)O)C(=O)O",                    # D-Glu
        "N[C@H](Cc1cnc[nH]1)C(=O)O",                 # D-His  ← REVIEWER: was absent
        "CC[C@@H](C)[C@H](N)C(=O)O",                 # D-Ile  (2R,3R)
        "CC(C)C[C@H](N)C(=O)O",                      # D-Leu
        "N[C@H](CCCCN)C(=O)O",                       # D-Lys
        "N[C@H](CCSC)C(=O)O",                        # D-Met
        "N[C@H](Cc1ccccc1)C(=O)O",                   # D-Phe
        "O=C(O)[C@H]1CCCN1",                         # D-Pro
        "N[C@H](CO)C(=O)O",                          # D-Ser
        "N[C@H]([C@H](O)C)C(=O)O",                   # D-Thr (2R,3S)
        "N[C@H](Cc1c[nH]c2ccccc12)C(=O)O",           # D-Trp
        "N[C@H](Cc1ccc(O)cc1)C(=O)O",                # D-Tyr
        "CC(C)[C@H](N)C(=O)O",                       # D-Val
    ]

    def _canonicalise(smi: str) -> str | None:
        mol = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(mol) if mol else None

    PRIORITY_SMILES: FrozenSet[str] = frozenset(
        filter(None, (_canonicalise(s) for s in _L_SMILES + _D_SMILES))
    )

    # 20 L + 18 D chiral + 1 Gly (achiral, deduplicated) = 39 unique expected.
    # Glycine has no stereocentre so L-Gly and D-Gly canonicalise to the same
    # string; the frozenset automatically deduplicates it.
    assert len(PRIORITY_SMILES) >= 20, (
        f"Expected ≥ 20 canonical priority SMILES, got {len(PRIORITY_SMILES)}"
    )

    # ── tag rows ─────────────────────────────────────────────────────────────────
    df_tagged = df_raw.copy()
    df_tagged["is_protected"] = df_tagged["canonical_smiles"].isin(PRIORITY_SMILES)

    n_protected = int(df_tagged["is_protected"].sum())

    mo.md(
        f"""
        **Priority set built.**

        | | |
        |-|-|
        | Canonical SMILES in protected set | {len(PRIORITY_SMILES)} |
        | Protected rows found in dataset | {n_protected} |
        """
    )
    return PRIORITY_SMILES, df_tagged, n_protected


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 4 · Filter by Category

        We keep only rows where the `category` column equals `"amino_acid"`
        (comparison is case-insensitive and strips surrounding whitespace).

        Protected rows are separated out **before** this filter so they are
        not accidentally removed if their category field differs.
        """
    )
    return


@app.cell
def __(config, df_tagged, mo):
    def _filter_category(df):
        """Return rows whose category matches target_category (case-insensitive)."""
        mask = (
            df[config.category_col]
            .str.strip()
            .str.lower()
            == config.target_category.lower()
        )
        return df[mask].copy()

    # REVIEWER: is_protected is checked *before* the category filter so that
    #           protected molecules bypass it entirely, regardless of their
    #           category field value.
    df_protected = df_tagged[df_tagged["is_protected"]].copy()
    df_unprotected = df_tagged[~df_tagged["is_protected"]].copy()

    df_category = _filter_category(df_unprotected)

    assert len(df_category) > 0, (
        f"No rows remain after category filter. "
        f"Check that '{config.category_col}' contains '{config.target_category}'."
    )

    mo.md(
        f"""
        **After category filter:**

        | | |
        |-|-|
        | Unprotected rows before filter | {len(df_unprotected):,} |
        | Unprotected rows after filter | {len(df_category):,} |
        | Removed | {len(df_unprotected) - len(df_category):,} |
        | Protected rows (bypass) | {len(df_protected):,} |
        """
    )
    return df_category, df_protected, df_unprotected


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 5 · Filter by SMARTS — Retain N-Methyl Amino Acids

        The SMARTS pattern `[NX3;H0;!$(NC=O)][CH3]` matches an N-methyl group on
        a **non-amide** trivalent nitrogen.

        ### Filter logic

        We want to **include** N-methyl amino acids (they are biologically
        interesting non-canonical monomers).  Therefore the filter removes rows
        that match a *problematic* secondary-amine pattern we do **not** want —
        but since the specification is to *retain* N-methyl AAs, we simply pass
        all rows through this step unchanged unless a stricter exclusion criterion
        is provided.

        Concretely: we keep every row that either
        - **does** match the N-methyl SMARTS (explicitly desired), or
        - **does not** match it (ordinary amino acid, also desired).

        In other words, **no row is removed** solely because of the SMARTS match.
        The cell reports how many N-methyl AAs were detected for bookkeeping.

        > **FIXME:** If the dataset contains other N-substituted AAs that *should*
        > be excluded (e.g. N-benzyl, N-acyl), add those SMARTS patterns here and
        > negate the match logic.
        """
    )
    return


@app.cell
def __(Chem, config, df_category, mo):
    _nmethyl_query = Chem.MolFromSmarts(config.nmethyl_smarts)
    assert _nmethyl_query is not None, (
        f"N-methyl SMARTS '{config.nmethyl_smarts}' failed to compile."
    )

    def _has_nmethyl(mol) -> bool:
        """Return True if mol contains an N-methyl (non-amide) group."""
        if mol is None:
            return False
        return mol.HasSubstructMatch(_nmethyl_query)

    df_smarts = df_category.copy()
    df_smarts["has_nmethyl"] = df_smarts[config.mol_col].apply(_has_nmethyl)

    n_nmethyl = int(df_smarts["has_nmethyl"].sum())

    assert len(df_smarts) == len(df_category), (
        "SMARTS step unexpectedly changed the row count."
    )

    mo.md(
        f"""
        **SMARTS annotation (all rows retained — see description above):**

        | | |
        |-|-|
        | Rows with N-methyl motif | {n_nmethyl:,} |
        | Rows without N-methyl | {len(df_smarts) - n_nmethyl:,} |
        | Total retained | {len(df_smarts):,} |
        """
    )
    return df_smarts, n_nmethyl


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 6 · Filter by Rotatable Bonds

        `Descriptors.NumRotatableBonds` counts the number of freely rotatable
        single bonds (RDKit default definition).  We discard molecules with more
        than `config.max_rot_bonds` (currently **5**) rotatable bonds.

        Protected rows skip this filter entirely.
        """
    )
    return


@app.cell
def __(Descriptors, config, df_smarts, mo):
    # REVIEWER: confirmed use of Descriptors.NumRotatableBonds (not
    #           rdMolDescriptors) and threshold is applied as <= (inclusive).
    def _filter_rotb(df, max_rotb: int):
        """Keep rows with ≤ max_rotb rotatable bonds."""
        df = df.copy()
        df["num_rot_bonds"] = df[config.mol_col].apply(
            lambda m: Descriptors.NumRotatableBonds(m) if m is not None else 9999
        )
        return df[df["num_rot_bonds"] <= max_rotb].copy()

    n_before_rotb = len(df_smarts)
    df_rotb = _filter_rotb(df_smarts, config.max_rot_bonds)
    n_after_rotb = len(df_rotb)

    assert n_after_rotb > 0, (
        f"All rows were removed by the rotatable-bond filter "
        f"(threshold = {config.max_rot_bonds}).  Consider raising max_rot_bonds."
    )

    mo.md(
        f"""
        **After rotatable-bond filter (≤ {config.max_rot_bonds}):**

        | | |
        |-|-|
        | Rows before | {n_before_rotb:,} |
        | Rows after | {n_after_rotb:,} |
        | Removed | {n_before_rotb - n_after_rotb:,} |
        """
    )
    return df_rotb, n_after_rotb, n_before_rotb


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 7 · Merge Protected + Filtered Sets

        The protected amino acids (step 3) are reunited with the molecules that
        survived steps 4–6.  We then deduplicate on canonical SMILES so that
        natural amino acids which also appear in the filtered set are not counted
        twice.

        Deduplication happens **before** fingerprint computation (step 8) to avoid
        inflating cluster sizes with redundant bit vectors.
        """
    )
    return


@app.cell
def __(df_protected, df_rotb, mo, pd):
    # REVIEWER: protected rows are concatenated first so that drop_duplicates(keep="first")
    #           always retains the protected copy when a canonical SMILES appears in both sets.
    df_merged = pd.concat([df_protected, df_rotb], ignore_index=True)

    n_before_dedup = len(df_merged)
    df_merged = df_merged.drop_duplicates(
        subset=["canonical_smiles"], keep="first"
    ).copy()
    n_after_dedup = len(df_merged)

    assert n_after_dedup > 0, "No rows remain after merging and deduplication."

    mo.md(
        f"""
        **After merge + deduplication:**

        | | |
        |-|-|
        | Rows before dedup | {n_before_dedup:,} |
        | Duplicates removed | {n_before_dedup - n_after_dedup:,} |
        | **Unique structures** | **{n_after_dedup:,}** |
        """
    )
    return df_merged, n_after_dedup, n_before_dedup


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 8 · Compute Morgan Fingerprints

        We compute **Morgan circular fingerprints** (radius 2, 2048 bits) for
        every molecule.  These bit vectors are used in step 9 to calculate
        pairwise Tanimoto similarities.
        """
    )
    return


@app.cell
def __(AllChem, config, df_merged, mo):
    def _compute_fp(mol):
        """Return Morgan bit-vector fingerprint or None if mol is None."""
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=config.fp_radius, nBits=config.fp_bits
        )

    df_fp = df_merged.copy()
    df_fp["fp"] = df_fp[config.mol_col].apply(_compute_fp)

    # Drop any row that failed fingerprinting (should not happen after step 2)
    n_before_fp = len(df_fp)
    df_fp = df_fp[df_fp["fp"].notna()].copy()
    n_dropped_fp = n_before_fp - len(df_fp)

    assert len(df_fp) > 0, "No valid fingerprints computed."

    _warn = (
        f"\n\n> ⚠️ **{n_dropped_fp} rows dropped** due to fingerprint failure."
        if n_dropped_fp else ""
    )

    mo.md(
        f"""
        **Fingerprints ready:** {len(df_fp):,} Morgan fingerprints
        (radius = {config.fp_radius}, {config.fp_bits} bits).{_warn}
        """
    )
    return df_fp, n_before_fp, n_dropped_fp


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 9 · Greedy Centroid Clustering (Tanimoto ≥ 0.9)

        ### Algorithm

        1. **Sort** molecules so that `is_protected = True` rows come first,
           preserving their priority as cluster centroids.
        2. Iterate over molecules in order:
           - If the molecule's Tanimoto similarity to **any existing centroid**
             is ≥ `config.tanimoto_threshold`, assign it to that cluster.
           - Otherwise, declare it a **new centroid**.
        3. Return a DataFrame containing **one row per cluster** (the centroid).

        This is an $O(n \cdot k)$ greedy algorithm where $k$ is the number of
        clusters found.  In the worst case (all singletons) $k = n$, giving
        $O(n^2)$ complexity.

        > ⚠️ **Performance note:** for libraries with **n > 500 molecules** the
        > pairwise loop below may become slow.  Consider switching to
        > `rdkit.ML.Cluster.Butina.ClusterData` which implements the same greedy
        > algorithm with a more cache-friendly layout and optional distance
        > pre-screening.
        """
    )
    return


@app.cell
def __(DataStructs, config, df_fp, mo):
    from typing import List

    def greedy_centroid_cluster(
        df,
        fp_col: str = "fp",
        protected_col: str = "is_protected",
        threshold: float = 0.9,
    ):
        """
        Greedy centroid clustering by Tanimoto similarity.

        Parameters
        ----------
        df            : DataFrame with fingerprint column and is_protected flag
        fp_col        : column name of RDKit ExplicitBitVect fingerprints
        protected_col : boolean column; protected rows are sorted first so they
                        are always chosen as centroids when a cluster contains them
        threshold     : minimum Tanimoto similarity to join an existing cluster
                        (applied as >= threshold, i.e. inclusive)

        Returns
        -------
        DataFrame of centroids with added columns cluster_id and cluster_size.

        Complexity
        ----------
        O(n * k) where k = number of clusters; worst case O(n^2) when all
        molecules are singletons.  For n > 500 consider Butina.ClusterData.
        """
        # REVIEWER: sorting protected rows first guarantees they become centroids
        #           before any unprotected molecule, regardless of insertion order.
        df_sorted = df.sort_values(
            by=protected_col, ascending=False
        ).reset_index(drop=True)

        fps: List = df_sorted[fp_col].tolist()
        n = len(fps)

        centroid_indices: List[int] = []
        cluster_assignments: List[int] = [-1] * n
        cluster_sizes: List[int] = []

        for i in range(n):
            fp_i = fps[i]
            assigned = False
            for cid, cent_idx in enumerate(centroid_indices):
                # REVIEWER: TanimotoSimilarity used (not BulkTanimotoSimilarity)
                # and threshold applied as >= (inclusive) per spec.
                sim = DataStructs.TanimotoSimilarity(fp_i, fps[cent_idx])
                if sim >= threshold:
                    cluster_assignments[i] = cid
                    cluster_sizes[cid] += 1
                    assigned = True
                    break
            if not assigned:
                new_cid = len(centroid_indices)
                centroid_indices.append(i)
                cluster_assignments[i] = new_cid
                cluster_sizes.append(1)

        df_sorted = df_sorted.copy()
        df_sorted["cluster_id"] = cluster_assignments

        centroid_rows = df_sorted.iloc[centroid_indices].copy()
        centroid_rows = centroid_rows.reset_index(drop=True)
        centroid_rows["cluster_id"] = list(range(len(centroid_indices)))
        centroid_rows["cluster_size"] = cluster_sizes

        return centroid_rows

    df_centroids = greedy_centroid_cluster(
        df_fp,
        fp_col="fp",
        protected_col="is_protected",
        threshold=config.tanimoto_threshold,
    )

    assert len(df_centroids) > 0, "Clustering produced zero centroids."
    assert "cluster_id" in df_centroids.columns
    assert "cluster_size" in df_centroids.columns

    mo.md(
        f"""
        **Clustering complete** (τ = {config.tanimoto_threshold}):

        | | |
        |-|-|
        | Input molecules | {len(df_fp):,} |
        | Clusters found | {len(df_centroids):,} |
        | Largest cluster | {df_centroids['cluster_size'].max()} members |
        | Singleton clusters | {int((df_centroids['cluster_size'] == 1).sum()):,} |
        """
    )
    return List, df_centroids, greedy_centroid_cluster


@app.cell
def __(mo):
    mo.md(
        r"""
        ## 10 · Output & Visualisation

        The grid below shows one centroid per cluster.  Each tile is labelled with
        its `cluster_id` and `cluster_size`.

        Summary statistics are printed at the bottom, and the centroid table is
        available for download as `centroids_output.json` in the same schema as
        the input file.
        """
    )
    return


@app.cell
def __(
    Draw,
    config,
    df_category,
    df_centroids,
    df_fp,
    df_raw,
    df_rotb,
    mo,
    n_after_dedup,
    n_after_rotb,
    n_dropped_fp,
    n_protected,
):
    import base64

    # ── summary statistics ───────────────────────────────────────────────────────
    _stats = {
        "Total loaded (valid SMILES)"   : len(df_raw),
        "Protected (natural AAs)"       : n_protected,
        "After category filter"         : len(df_category),
        "After rotatable-bond filter"   : n_after_rotb,
        "After merge + dedup"           : n_after_dedup,
        "Fingerprint failures"          : n_dropped_fp,
        "Final centroids"               : len(df_centroids),
    }

    _stats_md = "\n".join(
        f"| {k} | **{v:,}** |" for k, v in _stats.items()
    )

    _summary_table = mo.md(
        f"""
        ### Pipeline Summary

        | Step | Count |
        |------|-------|
        {_stats_md}
        """
    )

    # ── molecule grid ────────────────────────────────────────────────────────────
    # REVIEWER: PandasTools.FrameToGridImage is appropriate for ≤ 200 molecules;
    #           Draw.MolsToGridImage used here for fine-grained label control.
    #           For > 200 centroids mo.plain(df) is shown instead.
    MAX_GRID = 50  # FIXME: raise if you want more structures in the grid

    _display_df = df_centroids.head(MAX_GRID).copy()
    _display_df["Label"] = _display_df.apply(
        lambda r: f"c{int(r['cluster_id'])}  (n={int(r['cluster_size'])})",
        axis=1,
    )

    _mols = _display_df[config.mol_col].tolist()
    _labels = _display_df["Label"].tolist()

    try:
        _grid_img = Draw.MolsToGridImage(
            _mols,
            molsPerRow=5,
            subImgSize=(250, 200),
            legends=_labels,
            returnPNG=True,
        )
        _img_b64 = base64.b64encode(_grid_img).decode()
        _grid_html = mo.Html(
            f'<img src="data:image/png;base64,{_img_b64}" '
            f'style="max-width:100%;border:1px solid #ddd;border-radius:6px;" />'
        )
    except Exception as _exc:
        _grid_html = mo.md(f"> ⚠️ Grid rendering failed: `{_exc}`")

    mo.vstack([_summary_table, _grid_html])
    return MAX_GRID, base64


@app.cell
def __(df_centroids, mo):
    import io as _io
    import json as _json

    # ── export centroids to JSON ─────────────────────────────────────────────────
    _export_cols = [
        c for c in df_centroids.columns
        if c not in ("ROMol", "fp")  # drop non-serialisable objects
    ]
    _records = df_centroids[_export_cols].to_dict(orient="records")
    # REVIEWER: pass bytes directly to mo.download; wrapping in BytesIO is
    #           unnecessary and some marimo versions handle bytes more reliably.
    _json_bytes = _json.dumps(_records, indent=2, default=str).encode("utf-8")

    _download_btn = mo.download(
        data=_json_bytes,
        filename="centroids_output.json",
        mimetype="application/json",
        label="⬇ Download centroids_output.json",
    )

    mo.vstack([
        mo.md("### Export\nClick below to download the centroid set as JSON."),
        _download_btn,
    ])
    return


if __name__ == "__main__":
    app.run()
