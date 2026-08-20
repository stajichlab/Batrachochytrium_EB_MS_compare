# Environments & Installations

## Primary Environment

- **Manager**: pixi (v0.71.3, `/rhome/jstajich/.pixi/bin/pixi`)
- **Python version**: 3.12.*
- **Created**: 2026-08-19

### Setup from scratch

```bash
pixi install
```

### Run the analysis pipeline

```bash
pixi run build-ordination-table
pixi run pcoa-ordination
pixi run differential-features
```

## Dependencies

See `pixi.toml` / `pixi.lock` (tracked). Conda-forge + bioconda, linux-64:
python 3.12, pandas, numpy, scipy, matplotlib, seaborn, scikit-learn.

## System Dependencies

- `python3.12` (system; needed for mycelium's own scripts, e.g.
  `validate_structure.py`, `generate_index.py` — plain-text yaml fallbacks).
- `pixi` on PATH for the project env.
- SIRIUS 6.3.12 used by `analysis/sirius_annotation/` (native run), invoked
  via the SLURM array script; not installed into the pixi env.
