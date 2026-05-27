"""build_figures.py

Run every figure-generation script in ``code/`` and copy the resulting PDF
artifacts into ``manuscript/figures/`` so the LaTeX build can pick them up.

A "figure script" is any file in ``code/`` whose name matches the pattern
``calc-*.py``.  Each script is invoked as

    python <script> --out <code/figures/figure_X> --format pdf,png

and is expected to honour the ``--out`` and ``--format`` arguments.

Usage:
    python scripts/build_figures.py            # build every calc-*.py
    python scripts/build_figures.py calc-01    # build just calc-01.py
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CODE_DIR = REPO_ROOT / "code"
CODE_OUT = CODE_DIR / "figures"
MANUSCRIPT_FIGURES = REPO_ROOT / "manuscript" / "figures"


def discover_scripts(filters: list[str]) -> list[Path]:
    """Return the figure scripts to run, optionally filtered by stem."""
    all_scripts = sorted(CODE_DIR.glob("calc-*.py"))
    if not filters:
        return all_scripts
    selected: list[Path] = []
    for f in filters:
        match = next((s for s in all_scripts if s.stem == f or s.name == f), None)
        if match is None:
            raise SystemExit(f"error: no script matching '{f}' in {CODE_DIR}")
        selected.append(match)
    return selected


def figure_index(script: Path) -> str:
    """Extract a numeric figure index from a script name.

    ``calc-01.py``       -> ``1``
    ``calc-2-extra.py``  -> ``2`` (descriptive suffixes after the number are stripped)
    ``calc-foo.py``      -> ``foo`` (non-numeric stems are passed through verbatim)
    """
    stem = script.stem
    if not stem.startswith("calc-"):
        return stem
    tail = stem.split("-", 1)[1]
    head = tail.split("-", 1)[0]
    if head.isdigit():
        return str(int(head))
    return tail


def run_script(script: Path) -> Path:
    """Execute one figure script and return the produced PDF path."""
    idx = figure_index(script)
    out_basename = CODE_OUT / f"figure_{idx}"
    out_basename.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(script),
        "--out",
        str(out_basename),
        "--format",
        "pdf,png",
    ]
    print(f"[figures] running: {' '.join(cmd)}")
    subprocess.check_call(cmd, cwd=CODE_DIR)
    return out_basename.with_suffix(".pdf")


def copy_to_manuscript(pdf: Path) -> Path:
    MANUSCRIPT_FIGURES.mkdir(parents=True, exist_ok=True)
    target = MANUSCRIPT_FIGURES / pdf.name
    shutil.copy2(pdf, target)
    png = pdf.with_suffix(".png")
    if png.exists():
        shutil.copy2(png, MANUSCRIPT_FIGURES / png.name)
    print(f"[figures] copied -> {target.relative_to(REPO_ROOT)}")
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "scripts",
        nargs="*",
        help="Optional list of script stems (e.g. 'calc-01') to build. Defaults to all.",
    )
    args = parser.parse_args(argv)

    scripts = discover_scripts(args.scripts)
    if not scripts:
        print("[figures] no calc-*.py scripts found, nothing to do.")
        return 0

    for s in scripts:
        pdf = run_script(s)
        copy_to_manuscript(pdf)

    print(f"[figures] done -- {len(scripts)} figure(s) ready in {MANUSCRIPT_FIGURES.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
