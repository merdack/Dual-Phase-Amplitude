# Dual Phase-Amplitude Coupling in Synchronized Scalar Fields

> **A Dynamical Phase-Stress Approach** — manuscript for [*Physical Review D*](https://journals.aps.org/prd/) (`revtex4-2`, `prd`, two-column).

This repository holds the modular LaTeX source, reproducible Python figure scripts, and build pipeline for a covariant **phase-locked Klein–Gordon (PLKG)** model: $N$ complex scalars in flat FLRW spacetime with an $S_N$-symmetric democratic interaction whose polar form yields Kuramoto-type phase couplings, coupled amplitude–phase equations, the stress–energy tensor, stability analysis, and illustrative numerics. The layout supports parallel edits to sections, appendices, figures, or bibliography without merge conflicts.

**Output PDF:** `build/main.pdf` (not at the repository root).

---

## Table of contents

- [What is in the paper](#what-is-in-the-paper)
- [Project layout](#project-layout)
- [Quick start](#quick-start)
- [Building the manuscript](#building-the-manuscript)
- [Physical Review D conventions](#physical-review-d-conventions)
- [Optional draft overlays](#optional-draft-overlays)
- [Regenerating figures](#regenerating-figures)
- [How to contribute](#how-to-contribute)
- [License](#license)

---

## What is in the paper

| Part | File(s) | Content |
|------|---------|---------|
| Introduction | `sections/01_introduction.tex` | Motivation; PLKG vs standard KG and Kuramoto networks |
| Theory | `sections/02_theory.tex` | Lagrangian, EoM, polar split, $T_{\mu\nu}$, continuum limit |
| Results | `sections/03_results.tex` | Analytical bullets; Figs. 1–6 and 9 (homogeneous ODE, FRW schematic, stability, 3D grid) |
| Discussion | `sections/04_discussion.tex` | FDM, BBN/CMB scalings, structure formation, toy TOV (Figs. 7–8, 10) |
| Conclusion | `sections/05_conclusion.tex` | Summary |
| Appendix A | `appendices/A_conservation.tex` | $\nabla_\mu T^{\mu\nu}=0$ on shell |
| Appendix B | `appendices/B_stability.tex` | Linearization, Jeans scale, $K_{\mathrm{cosmo}}\sim 3H^2N$ |
| Appendix C | `appendices/C_reproducibility.tex` | How to rerun `calc-*.py` and rebuild the PDF |

Title, author, and abstract live in [`manuscript/metadata.tex`](manuscript/metadata.tex). Prose is written in **passive voice**; equations and section structure are unchanged from the modular split.

---

## Project layout

```
dual-phase-amplitude-01/
├── README.md                 <- you are here
├── AGENTS.md                 <- instructions for AI coding agents
├── CONTRIBUTING.md           <- branch workflow, LaTeX/citation rules
├── LICENSE
├── Makefile                  <- make | make figures | make clean
├── scripts/
│   ├── build.ps1             <- Windows / PowerShell build
│   ├── build.sh              <- POSIX bash build
│   └── build_figures.py      <- runs code/calc-*.py -> manuscript/figures/
├── build/                    <- LaTeX output (main.pdf, aux, .bbl)
├── manuscript/
│   ├── main.tex              <- document class + \\input chain only
│   ├── metadata.tex          <- title, author, abstract
│   ├── preamble/
│   │   ├── packages.tex      <- packages + \\bibliographystyle{apsrev4-2}
│   │   └── macros.tex        <- PLKG macros; \\Refcite / \\Refscite; \\showdraft
│   ├── sections/             <- 01_introduction .. 05_conclusion
│   ├── drafts/               <- optional WIP (\\showdrafttrue only)
│   ├── appendices/           <- A_conservation, B_stability, C_reproducibility
│   ├── bibliography/
│   │   └── references.bib    <- BibTeX (apsrev4-2)
│   └── figures/              <- figure_1.pdf .. figure_10.pdf (from build_figures.py)
├── code/
│   ├── README.md             <- figure catalogue (calc-01 .. calc-10)
│   ├── requirements.txt
│   ├── calc-01.py .. calc-10.py
│   └── figures/              <- generated pdf/png (often git-ignored)
├── docs/                     <- notes not compiled into the PDF
└── .github/workflows/        <- CI: figures + PDF on push
```

The article root is **`manuscript/main.tex`**. An older monolithic `main.tex` at the repo root, if present locally, is obsolete.

---

## Quick start

**Prerequisites**

- TeX: `pdflatex` and `bibtex` (TeX Live, MiKTeX, etc.)
- Python 3.10+

```bash
git clone <repo-url> dual-phase-amplitude-01
cd dual-phase-amplitude-01

python -m pip install -r code/requirements.txt
make                              # POSIX
# or
pwsh scripts/build.ps1            # Windows / PowerShell 7+
```

Open **`build/main.pdf`**.

Text-only iteration (no figure regeneration):

```bash
bash scripts/build.sh --skip-figures
pwsh scripts/build.ps1 -SkipFigures
```

---

## Building the manuscript

The pipeline runs in order:

1. **Figures** — `scripts/build_figures.py` executes each `code/calc-*.py` with `--out code/figures/figure_<N>` and copies `figure_<N>.pdf` into `manuscript/figures/`.
2. **LaTeX** — `pdflatex` → `bibtex` → `pdflatex` ×2, with `-output-directory=build` and source `manuscript/main.tex`.
3. **Cleanup (optional)** — drop aux/log files while keeping the PDF.

| Goal | POSIX | PowerShell |
|------|--------|------------|
| Full build | `make` | `pwsh scripts\build.ps1` |
| Rebuild figures + PDF | `bash scripts/build.sh --figures` | `pwsh scripts\build.ps1 -Figures` |
| PDF only | `bash scripts/build.sh --skip-figures` | `pwsh scripts\build.ps1 -SkipFigures` |
| Figures only | `make figures` | `python scripts/build_figures.py` |
| Clean aux, keep PDF | `make clean` | `pwsh scripts\build.ps1 -Clean -KeepPdf` |

Manual LaTeX (from repo root):

```bash
latexmk -pdf -outdir=build manuscript/main.tex
```

Before merging, confirm the log has no undefined citations or multiply-defined labels.

---

## Physical Review D conventions

The manuscript targets APS submission style:

| Item | Where |
|------|--------|
| Document class | `revtex4-2` with `aps`, `prd`, `twocolumn`, `superscriptaddress`, `nofootinbib`, `longbibliography` |
| BibTeX style | `apsrev4-2` in `manuscript/preamble/packages.tex` (**before** `\begin{document}`) |
| Database | `manuscript/bibliography/references.bib` — keys `LastNameYear`, APS journal abbreviations, `doi` when available |
| In-text cites | Non-breaking space: `model~\cite{Guth1981}` |
| Inline “Ref. 5” | `Ref.~\onlinecite{key}` or macros `\Refcite{key}` / `\Refscite{key1,key2}` in `macros.tex` |

Full citation rules are in [`CONTRIBUTING.md`](CONTRIBUTING.md) §3. Do not put `%` comment lines inside `.bib` files (BibTeX treats `%` as data); use `@comment{ ... }` instead.

---

## Optional draft overlays

Work-in-progress blocks under `manuscript/drafts/` compile only when draft mode is on:

1. In `manuscript/main.tex`, uncomment `\showdrafttrue` immediately after `\input{preamble/macros}` (default is off; CI stays clean).
2. Example overlay: `drafts/02_theory_sketch.tex` is inserted after `sections/02_theory` inside `\ifshowdraft ... \fi`.
3. Use unique `\label{...}` keys in draft files so they never clash with the main text.

Long notes that are never part of the PDF belong under [`docs/`](docs/) (e.g. `docs/physics_notes.md`).

---

## Regenerating figures

Ten scripts, `calc-01.py` through `calc-10.py`, produce `figure_1` … `figure_10` for the Results and Discussion sections. Each script:

- accepts `--out <basename>` and writes `.pdf` / `.png`;
- accepts `--prd-width column|full` for PRD float widths;
- seeds RNGs for reproducible output.

```bash
python scripts/build_figures.py           # all figures
python scripts/build_figures.py calc-03   # one figure
python code/calc-01.py --show             # interactive preview
```

See [`code/README.md`](code/README.md) for the script-to-section mapping and rules for new figures.

---

## How to contribute

Read [`CONTRIBUTING.md`](CONTRIBUTING.md) for branch prefixes (`sec/`, `app/`, `fig/`, `bib/`, `meta/`, `build/`), commit style, and the PR checklist. Automated agents should follow [`AGENTS.md`](AGENTS.md).

**TL;DR**

- One logical change per branch (one section, one appendix, one figure, or one bib entry).
- Edit prose in `manuscript/sections/*.tex` and `manuscript/appendices/*.tex`, not in `main.tex` (except new `\input` lines or `\showdrafttrue`).
- Run `make` or `build.ps1` before opening a PR; keep CI green.

---

## License

Python code under `code/` and `scripts/` is [MIT](LICENSE). Manuscript text and figures are copyright the author(s); public release terms are CC BY 4.0 as stated in the repository license notes.
