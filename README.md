# dutils0

Utility package.

## Installation

### pip (editable)

```bash
pip install -e .
```

### Pixi

[Pixi](https://pixi.sh) uses this repo’s `pyproject.toml` and installs the package in editable mode. From the repo root:

```bash
pixi install
```

Then run Python in the pixi environment:

```bash
pixi run python
```

Or activate the environment and use your normal workflow:

```bash
pixi shell
python -c "import dutils0; print(dutils0.__version__)"
```
