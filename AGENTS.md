# AGENTS.md — AI agent instructions

This file tells automated coding agents (Cursor, Copilot, etc.) how to work
safely and productively in **dual-phase-amplitude-01**: a LaTeX manuscript
for *Physical Review D* (revtex4-2), plus reproducible Python figures and a
small build pipeline.

Human-oriented workflow details live in [`CONTRIBUTING.md`](CONTRIBUTING.md)
and [`README.md`](README.md). Follow those for branches and PR scope; use
this file for **what to touch, how to build, and what to avoid**.

---

## 1. Project purpose

- **Manuscript:** Klein–Gordon fields with Kuramoto-like phase coupling;
  amplitude–phase (polar) equations; stress–energy tensor; appendices on
  conservation and stability.
- **Stack:** `pdflatex` + `bibtex`, `revtex4-2` (`prd`, two-column), Python
  3.10+ (NumPy / SciPy / Matplotlib) for figures.
- **Output PDF:** `build/main.pdf` (LaTeX writes to `build/`; do not assume
  PDF at repo root).

---

## 2. Repository layout (edit here)

| Path | Purpose |
|------|---------|
| `manuscript/main.tex` | Master file: `\input` only. **Do not put prose here.** |
| `manuscript/metadata.tex` | Title, authors, affiliations, abstract |
| `manuscript/preamble/packages.tex` | Package imports |
| `manuscript/preamble/macros.tex` | Shared `\providecommand` macros only |
| `manuscript/sections/0X_*.tex` | Main-body sections (introduction → conclusion) |
| `manuscript/appendices/*.tex` | Appendices (single `\appendix` lives in `main.tex`) |
| `manuscript/bibliography/references.bib` | BibTeX database (`apsrev4-2` style) |
| `manuscript/figures/*.pdf` | Figures consumed by LaTeX (often copied from build) |
| `code/calc-*.py` | One script per figure; reproducible, CLI-driven |
| `code/requirements.txt` | Python dependencies |
| `scripts/build.ps1`, `scripts/build.sh` | Full PDF build |
| `scripts/build_figures.py` | Runs all `code/calc-*.py`, copies PDFs into `manuscript/figures/` |
| `Makefile` | POSIX entry point (`make`, `make figures`, `make clean`) |
| `docs/physics_notes.md` | Extra notes / derivations (not part of the compiled PDF unless linked) |
| `.github/workflows/` | CI builds the PDF on push |

---

## 3. Golden rules for agents

1. **Scope:** Prefer a *minimal* diff: one logical change (one section, one
   appendix, one figure script, or one bibliography block). No drive-by
   refactors across unrelated files unless the user explicitly asks.

2. **Where prose and equations live:** Edit `manuscript/sections/*.tex`,
   `manuscript/appendices/*.tex`, or `manuscript/metadata.tex`. Touch
   `manuscript/main.tex` only to add/remove/reorder `\input{...}` lines.

3. **LaTeX style:** One sentence per line in `.tex` files. Use semantic
   `\label{eq:...}` / `\label{fig:...}`; avoid hard-coded equation numbers in
   prose — use `\eqref`.

4. **Macros:** New commands go in `manuscript/preamble/macros.tex` with
   `\providecommand` so they do not clash with revtex internals.

5. **Bibliography:** Add entries only to `manuscript/bibliography/references.bib`,
 sorted by cite key (`LastNameYear`). **BibTeX does not treat `%` as a
 comment inside `.bib` files** — use `@comment{ ... }` for notes, or avoid
 `%` lines entirely in `.bib`. Use APS abbreviations (`Phys. Rev. D`, etc.),
 include `doi` when available, and keep `\bibliographystyle{apsrev4-2}` in
 `preamble/packages.tex` (not after `\begin{document}`). Cite in prose as
 `~\cite{key}`; use `\onlinecite` / `\Refcite` / `\Refscite` when the number
 must appear inline per PRD convention.

6. **Language:** The manuscript is **English** (APS style). Keep terminology
   consistent with existing sections (e.g. KKG, FLRW, \(T_{\mu\nu}\)).

7. **Consistency:** After changing signs in \(\mathcal{L}\) or \(T_{\mu\nu}\),
   verify equations of motion and appendix derivations still agree. Run a
   full build and fix undefined references.

---

## 4. How to build

### Prerequisites

- TeX: `pdflatex` and `bibtex` (TeX Live / MiKTeX).
- Python 3.10+: `python -m pip install -r code/requirements.txt`.

### Commands

| Goal | POSIX | Windows (PowerShell) |
|------|--------|----------------------|
| Full build (figures + PDF) | `make` | `pwsh scripts\build.ps1` |
| Rebuild figures then PDF | `bash scripts/build.sh --figures` | `pwsh scripts\build.ps1 -Figures` |
| PDF only (skip figures) | `bash scripts/build.sh --skip-figures` | `pwsh scripts\build.ps1 -SkipFigures` |
| Figures only | `make figures` or `python scripts/build_figures.py` | `python scripts\build_figures.py` |
| Clean aux files, keep PDF | `make clean` (if defined) / see README | `pwsh scripts\build.ps1 -Clean -KeepPdf` |

Alternative manual LaTeX pass (from repo root):

```bash
latexmk -pdf -outdir=build manuscript/main.tex
```

After substantive edits, **run a full build** and ensure the log has no
undefined citations or multiply-defined labels.

---

## 5. Figures (`code/calc-*.py`)

When adding or changing a figure script:

1. Name: `calc-<NN>-<short-name>.py` (existing: e.g. `calc-01.py`).
2. CLI: support `--out <basename>` → write `<basename>.pdf` and
   `<basename>.png`; optional `--show` for interactive preview.
3. **Seed all RNGs** for reproducible outputs across machines.
4. Top-of-file **docstring** describing what the figure shows.
5. Register the script in the catalogue table in `code/README.md`.
6. Run `python scripts/build_figures.py` (or full `make` / `build.ps1`) so
   `manuscript/figures/` contains the PDF referenced by `\includegraphics`.

Python style: PEP 8, ~100 columns, type hints in new code; prefer explicit
matplotlib keyword arguments for readability.

---

## 6. Git / PR expectations (agents assisting humans)

- **Branch prefixes** (see `CONTRIBUTING.md`): `sec/`, `app/`, `fig/`, `bib/`,
  `meta/`, `build/`.
- **Commits:** imperative mood, ~70 character subject
  (e.g. `theory: clarify U(1) breaking in Kuramoto term`).
- **Before merge:** figures pipeline OK, `build/main.pdf` builds, CI green
  (`.github/workflows/build.yml`).

---

## 7. What agents should not do

- Do not replace the modular manuscript with a single giant `.tex` unless the
  user requests a deliberate migration.
- Do not commit generated junk unrelated to the paper (local editor backup
  files, personal `notes.tex` at root) unless the project already tracks that
  pattern.
- Do not invent numerical bounds (e.g. constraints on coupling constants)
  without derivation or citation — use placeholders clearly or omit.
- Do not strip `%` comments that carry physical or sign-convention warnings
  unless updating the surrounding math for consistency.

---

## 8. Quick reference: user request → likely files

| User goal | Primary files |
|-----------|----------------|
| Introduction / motivation | `manuscript/sections/01_introduction.tex` |
| Lagrangian, EoM, polar split, \(T_{\mu\nu}\) | `manuscript/sections/02_theory.tex` |
| Results, figures, dispersion bullets | `manuscript/sections/03_results.tex` |
| Cosmology / observations / limitations | `manuscript/sections/04_discussion.tex` |
| Closing summary | `manuscript/sections/05_conclusion.tex` |
| \(\nabla_\mu T^{\mu\nu}\) proof | `manuscript/appendices/A_conservation.tex` |
| Stability / linearization | `manuscript/appendices/B_stability.tex` |
| Title / abstract | `manuscript/metadata.tex` |
| New macro | `manuscript/preamble/macros.tex` |
| New citation | `manuscript/bibliography/references.bib` + first `\cite` in prose |
| Sync figure | `code/calc-*.py`, `manuscript/sections/03_results.tex` (caption path) |
| CI / build | `scripts/*`, `Makefile`, `.github/workflows/*` |

---

## 9. License note (for context only)

`code/` and `scripts/` are MIT-licensed; manuscript text and figures follow
author copyright / CC BY 4.0 terms as stated in `README.md`. Do not change
license files unless the repository owner asks.

---

When in doubt, read [`CONTRIBUTING.md`](CONTRIBUTING.md) and match existing
tone, notation, and structure in the nearest section before editing.
