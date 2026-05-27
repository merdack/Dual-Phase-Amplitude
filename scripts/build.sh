#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# build.sh -- Build the manuscript on POSIX (Linux, macOS, WSL).
# Compiles from manuscript/ so \input and \bibliography paths resolve; PDF
# lands in ../build/main.pdf.
# -----------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANUSCRIPT_DIR="${REPO_ROOT}/manuscript"
BUILD_DIR="${REPO_ROOT}/build"
JOBNAME="main"
MAIN_TEX="main.tex"
ENGINE="${ENGINE:-pdflatex}"
FORCE_FIGURES=0
SKIP_FIGURES=0
DO_CLEAN=0
KEEP_PDF=0

usage() { sed -n '2,12p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --figures)        FORCE_FIGURES=1 ;;
        --skip-figures)   SKIP_FIGURES=1 ;;
        --clean)          DO_CLEAN=1 ;;
        --keep-pdf)       KEEP_PDF=1 ;;
        --engine)         shift; ENGINE="$1" ;;
        -h|--help)        usage 0 ;;
        *) echo "unknown option: $1" >&2; usage 1 ;;
    esac
    shift
done

step() { printf "\n==> %s\n" "$1"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "error: '$1' is required" >&2; exit 1; }; }

need "${ENGINE}"
need bibtex

# --- figures -----------------------------------------------------------------
if [[ ${SKIP_FIGURES} -eq 0 ]]; then
    fig_dir="${MANUSCRIPT_DIR}/figures"
    if [[ ${FORCE_FIGURES} -eq 1 ]] || ! compgen -G "${fig_dir}/*.pdf" > /dev/null 2>&1; then
        step "Generating figures via build_figures.py"
        PY="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
        if [[ -z "${PY}" ]]; then
            echo "error: python3 or python is required for figures" >&2
            exit 1
        fi
        "${PY}" "${REPO_ROOT}/scripts/build_figures.py"
    else
        step "Figures already present (use --figures to force a rebuild)"
    fi
fi

# --- LaTeX -------------------------------------------------------------------
mkdir -p "${BUILD_DIR}"
BUILD_DIR_ABS="$(cd "${BUILD_DIR}" && pwd)"

export BIBINPUTS="${MANUSCRIPT_DIR}${BIBINPUTS:+:${BIBINPUTS}}"

run_engine() {
    ( cd "${MANUSCRIPT_DIR}" && "${ENGINE}" \
        -interaction=nonstopmode \
        -halt-on-error \
        -output-directory="${BUILD_DIR_ABS}" \
        "${MAIN_TEX}" )
}

step "Running ${ENGINE} (pass 1)"
run_engine

step "Running bibtex"
( cd "${BUILD_DIR_ABS}" && bibtex "${JOBNAME}" ) || \
    echo "  bibtex returned non-zero; check ${BUILD_DIR_ABS}/${JOBNAME}.blg"

step "Running ${ENGINE} (pass 2)"
run_engine

step "Running ${ENGINE} (pass 3, final cross-references)"
run_engine

PDF="${BUILD_DIR}/${JOBNAME}.pdf"
if [[ ! -f "${PDF}" ]]; then
    echo "error: build finished but ${PDF} was not produced." >&2
    exit 1
fi

printf "\nPDF built successfully: %s\n" "${PDF}"

if [[ ${DO_CLEAN} -eq 1 ]]; then
    if [[ ${KEEP_PDF} -eq 1 ]]; then
        find "${BUILD_DIR}" -mindepth 1 ! -name "${JOBNAME}.pdf" -exec rm -rf {} +
        echo "Cleaned intermediate files in ${BUILD_DIR} (kept ${JOBNAME}.pdf)"
    else
        rm -rf "${BUILD_DIR}"
        echo "Removed ${BUILD_DIR}"
    fi
fi
