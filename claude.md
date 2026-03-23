
🧪 claude.md — Cheminformatics Utility (RDKit +
Analysis)
Project Overview
This project is a Python-based cheminformatics toolkit built on RDKit for molecular
analysis, transformation, and data exploration.
The goal is to:
Build reusable utilities for molecular descriptors, fingerprints, and structure handling
Enable rapid exploratory data analysis of chemical datasets
Support scriptable workflows (CLI + notebooks)
Core Stack
Python 3.10+
RDKit (primary chemistry toolkit)
NumPy / pandas (data handling)
matplotlib / seaborn / plotly (visualization)
Optional: scikit-learn (ML), tqdm (progress), typer/click (CLI)
Key Concepts
ChatGPT Get Plus
3/23/26, 6:05 PM Cheminformatics Utility Guide
https://chatgpt.com/c/69c17266-4ca0-832d-9b36-d829d11c7bf3 1/5
Molecules are primarily handled as:
SMILES strings
RDKit Mol objects
Data is typically stored in:
pandas DataFrames (rows = molecules)
Expensive operations (e.g., conformers, fingerprints) should be cached where possible
use the code in rdsl folder for selection of rdkit molecules
Coding Guidelines
General
Prefer small, composable functions
Avoid monolithic scripts — split into modules:
io.py (loading/saving)
descriptors.py
fingerprints.py
filters.py
visualization.py
RDKit Best Practices
Always sanitize molecules unless explicitly skipped
Handle invalid SMILES gracefully (never crash pipelines)
Use canonical SMILES when deduplicating
Avoid recomputing properties unnecessarily
Example pattern:
Data Handling Conventions
Input datasets usually contain:
smiles column (required)
optional metadata columns
Add derived columns instead of overwriting:
def smiles_to_mol(smiles: str) -> Mol | None:
 mol = Chem.MolFromSmiles(smiles)
 if mol is None:
 return None
 return mol
Python Run
3/23/26, 6:05 PM Cheminformatics Utility Guide
https://chatgpt.com/c/69c17266-4ca0-832d-9b36-d829d11c7bf3 2/5
mol
mw
logp
fp
Common Tasks (AI Should Help With)
1. Descriptor Calculation
Molecular weight
LogP
TPSA
H-bond donors/acceptors
2. Fingerprints
Morgan (ECFP)
Bit vectors or numpy arrays
Similarity calculations (Tanimoto)
3. Filtering
Rule-of-five filters
Substructure matching (SMARTS)
Removing duplicates
4. Visualization
Distribution plots of descriptors
Chemical space plots (PCA/UMAP)
Molecule grids (RDKit drawing)
Performance Considerations
Use vectorized pandas operations where possible
Cache fingerprints and descriptors
Avoid repeated MolFromSmiles calls
Consider multiprocessing for large datasets
Error Handling
Never fail on a single bad molecule
3/23/26, 6:05 PM Cheminformatics Utility Guide
https://chatgpt.com/c/69c17266-4ca0-832d-9b36-d829d11c7bf3 3/5
Log or collect invalid SMILES
Functions should return None or NaN instead of crashing
Testing
Include small test datasets (10–100 molecules)
Validate:
SMILES parsing
Descriptor correctness
Fingerprint consistency
Style Preferences
Type hints are required
Use clear variable names ( mol , smiles , df )
Prefer readability over cleverness
Include docstrings for all public functions
Example Workflow
What the AI Should Do
Suggest idiomatic RDKit solutions
Help refactor scripts into reusable modules
Optimize slow cheminformatics workflows
Propose useful descriptors or analyses
Avoid unnecessary dependencies
df = load_smiles_csv("data.csv")
df["mol"] = df["smiles"].apply(smiles_to_mol)
df = df.dropna(subset=["mol"])
df["mw"] = df["mol"].apply(calc_mw)
df["fp"] = df["mol"].apply(morgan_fp)
plot_descriptor_distribution(df, "mw")
Python Run
3/23/26, 6:05 PM Cheminformatics Utility Guide
https://chatgpt.com/c/69c17266-4ca0-832d-9b36-d829d11c7bf3 4/5
Do you like this personality?
What the AI Should NOT Do
Do not introduce heavy frameworks unless justified
Do not assume GPU/Deep Learning unless requested
Do not ignore chemical validity issues
Do not silently drop data without explanation
Nice-to-Have Extensions
CLI interface for batch processing
Molecule standardization pipeline
Scaffold analysis
Similarity search utilities
Integration with SDF/SMILES/Parquet formats
If you want, I can tailor this further for a specific use case like:
drug discovery / ADMET filtering
reaction informatics
library enumeration
ML-ready feature pipelines
Those benefit from slightly different claude.md emphasis.
3/23/26, 6:05 PM Cheminformatics Utility Guide
https://chatgpt.com/c/69c17266-4ca0-832d-9b36-d829d11c7bf3 5/5