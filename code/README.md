# Numerics & figure generation

All numerical experiments and figure-generation scripts that back the
manuscript live in this directory. The top-level build script copies the
generated figures into `../manuscript/figures/`, where `\includegraphics`
picks them up.

## Layout

```
code/
├── README.md              <- this file
├── requirements.txt       <- pinned dependencies (NumPy / SciPy / Matplotlib)
├── calc-01.py             <- Figure 1: KKG point-dynamics synchronization
├── calc-02.py             <- Figure 2: order parameter vs time for several K
├── calc-03.py             <- Figure 3: R vs K scan and crossing-time curve
├── calc-04.py             <- Figure 4: cosmological $K$ vs $3H(z)^2N$ schematic
├── calc-05.py             <- Figure 5: m_eff^2 vs K and Jeans boundary map
├── calc-06.py             <- Figure 6: stable $\langle\rho\rangle$ vs $K$, $m_{\mathrm{eff}}^2$ slice, $\lambda$ regularization
├── calc-07.py             <- Figure 7: 3D CDM vs KKG density surfaces (schematic)
├── calc-08.py             <- Figure 8: FDM vs KKG |Phi|^2 toy + Var(max rho)
├── calc-09.py             <- Figure 9: 3D grid Eqs. (14)+(16); optional Poisson / external Phi
├── calc-10.py             <- Figure 10: toy TOV M-R; baseline EoS vs K-stiffened sync pressure
└── figures/               <- generated artifacts (pdf + png), git-ignored
```

## Quick start

```bash
# 1. Create an isolated environment (optional but recommended)
python -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\Activate.ps1

# 2. Install the dependencies
pip install -r requirements.txt

# 3. Generate every figure used by the manuscript
python ../scripts/build_figures.py     # invokes every calc-*.py in this folder
```

To regenerate a single figure interactively:

```bash
python calc-01.py --show
```

## Conventions for new figure scripts

1. Place new scripts under `code/` with a descriptive prefix, e.g.
   `calc-02-power-spectrum.py`.
2. Accept `--out <basename>` and `--format pdf,png` (see `calc-01.py`); the
   top-level `build_figures.py` driver passes these automatically.
   Figure scripts should also accept `--prd-width column|full` for
   Physical Review single-column vs two-column figure widths (inches).
3. Always write a deterministic seed when randomness is involved; the
   build pipeline assumes regenerated figures are bit-identical.
4. Default output formats should be **both** `pdf` (for the manuscript)
   and `png` (for previews / GitHub).
5. Add a one-line summary of the new figure to this README.

## Figure catalogue

| Script        | Output filename       | Used in                       |
|---------------|-----------------------|-------------------------------|
| `calc-01.py`  | `figures/figure_1.*`  | Sec. III (Results), Fig. 1    |
| `calc-02.py`  | `figures/figure_2.*`  | Sec. III (Results), Fig. 2    |
| `calc-03.py`  | `figures/figure_3.*`  | Sec. III (Results), Fig. 3    |
| `calc-04.py`  | `figures/figure_4.*`  | Sec. III (Results), Fig. 4    |
| `calc-05.py`  | `figures/figure_5.*`  | Sec. III (Results), Fig. 5    |
| `calc-06.py`  | `figures/figure_6.*`  | Sec. III (Results), Fig. 6    |
| `calc-07.py`  | `figures/figure_7.*`  | Sec. IV (Discussion), Fig. 7     |
| `calc-08.py`  | `figures/figure_8.*`  | Sec. IV (Discussion), Fig. 8      |
| `calc-09.py`  | `figures/figure_9.*`  | Sec. III (Results), Fig.~9; 3D grid Eqs.~(14)--(16) + optional gravity |
| `calc-10.py`  | `figures/figure_10.*` | Sec. IV (Discussion), Fig.~10; toy TOV mass--radius vs $K$ |
