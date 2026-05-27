# Contributing

Thanks for your interest in helping with this paper! The repository is
structured so that several authors can work on different parts of the
manuscript or the numerics without merge conflicts. The guidelines below
keep that in working order.

## 1. Pick a small unit of work

Choose **one** of the following granularities per branch / pull request:

- A single section (`manuscript/sections/0X_*.tex`)
- A single appendix (`manuscript/appendices/*.tex`)
- A single figure (one `code/calc-*.py`)
- A single bibliography entry
- A build / tooling improvement (scripts, CI, docs)

Avoid mega-PRs that touch many unrelated files.

## 2. Branching and commits

```bash
git checkout -b sec/results-numerical-illustration
# ... edit one or two files ...
git commit -m "results: tighten wording around fig. 1"
git push -u origin sec/results-numerical-illustration
```

Branch-name prefixes used in this repo:

| Prefix     | Use for                                            |
| ---------- | -------------------------------------------------- |
| `sec/`     | Edits to `manuscript/sections/*`                   |
| `app/`     | Edits to `manuscript/appendices/*`                 |
| `fig/`     | New figure scripts or visual changes to a figure   |
| `bib/`     | Additions to `manuscript/bibliography/references.bib` |
| `meta/`    | Title block, authors, abstract                     |
| `build/`   | Tooling: `scripts/`, `Makefile`, CI, `.gitignore`  |

Commits should be in the imperative present tense and stay under ~70 cols
for the subject (`results: rewrite predictions list`).

## 3. LaTeX conventions

- **Do not edit `manuscript/main.tex`** unless you are adding a new
  section or moving a structural piece. All prose belongs in the included
  files.
- One sentence per line in `.tex` files. This dramatically improves
  diff and review quality.
- Use semantic equation labels (`\label{eq:dispersion}`,
  `\label{eq:lambda_sync}`) rather than positional ones.
- New macros go in `manuscript/preamble/macros.tex` and must be
  `\providecommand`-ed so they cannot clash with revtex internals.
- Bibliographic entries go in `manuscript/bibliography/references.bib`,
  sorted alphabetically by cite key (`LastNameYear`). **BibTeX does not treat
  `%` as a comment inside `.bib` files** — use `@comment{ ... }` blocks for
  prose notes, or keep the file free of `%` lines.
- **Physical Review D citations:** the manuscript uses `revtex4-2` with
  `apsrev4-2` (set in `preamble/packages.tex` before `\begin{document}`) and
  the `longbibliography` class option so cited article titles appear in the
  reference list. In prose, attach citations with a non-breaking space:
  `...\ model~\cite{Guth1981}`. When the reference number must read inline
  (e.g., “Ref. 5”), use `Ref.~\onlinecite{key}` or `Refs.~\onlinecite{key1,key2}`
  (macros `\Refcite` / `\Refscite` in `macros.tex`). Combine multiple keys in
  one `\cite{key1,key2}`; use `\cite{key1,*key2}` only when BibTeX should merge
  distinct entries into a single bibliography item.

## 4. Figure conventions

Every figure script under `code/` must:

1. Sit in a single file named `calc-<NN>-<short-name>.py`.
2. Accept `--out <basename>` and write `<basename>.pdf` and
   `<basename>.png`.
3. Accept `--show` to optionally pop up an interactive window.
4. Seed every random source for byte-identical reproducibility.
5. Have a short docstring at the top describing the figure.
6. Be added to the catalogue table in `code/README.md`.

The build driver (`scripts/build_figures.py`) discovers any script
matching `code/calc-*.py` automatically, runs it, and copies the PDF into
`manuscript/figures/`.

## 5. Before opening a pull request

1. Run the figures: `python scripts/build_figures.py`.
2. Build the PDF: `make` (POSIX) or `pwsh scripts\build.ps1` (Windows).
3. Make sure `build/main.pdf` is up to date and free of unresolved
   references / undefined citations (check the warning summary at the
   end of the log).
4. If you added new Python packages, update
   `code/requirements.txt` with a `>=`-bounded entry.
5. CI must pass on your PR (see `.github/workflows/build.yml`).

## 6. Code style (Python)

- Target Python 3.10+; use type hints in new code.
- Follow [PEP 8](https://peps.python.org/pep-0008/) (4-space indent,
  ~100-col lines).
- Prefer explicit named arguments when calling matplotlib helpers --
  the figures must be readable for someone reviewing the script later.
- Keep numerical kernels in pure NumPy / SciPy whenever possible.

## 7. Reporting problems

If you find a derivation error, a typo, or a build problem, open an issue
that contains:

- A clear description of the expected and actual behaviour.
- The relevant file path and (for LaTeX) the line number.
- If a build problem: the failing command and the final ~30 lines of the
  `.log` file.

Thanks again for contributing!
