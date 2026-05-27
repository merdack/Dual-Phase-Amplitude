# =============================================================================
# Makefile -- thin wrapper over scripts/build.sh for convenience.
#
# Default target builds the PDF (rebuilding figures only when missing).
# All real logic lives in scripts/build.sh and scripts/build.ps1 -- this file
# is just a tab-and-comma friendly entry point.
# =============================================================================

REPO_ROOT     := $(CURDIR)
MANUSCRIPT    := $(REPO_ROOT)/manuscript
BUILD_DIR     := $(REPO_ROOT)/build
JOBNAME       := main
ENGINE        ?= pdflatex
PYTHON        ?= python3

PDF           := $(BUILD_DIR)/$(JOBNAME).pdf

.PHONY: all paper figures clean distclean help

all: paper

help:
	@echo "Available targets:"
	@echo "  make            -- build the PDF (rebuilds figures if missing)"
	@echo "  make figures    -- regenerate every figure in manuscript/figures"
	@echo "  make paper      -- compile manuscript/main.tex -> build/main.pdf"
	@echo "  make clean      -- remove intermediate LaTeX artifacts"
	@echo "  make distclean  -- remove build/ and generated figures entirely"

figures:
	$(PYTHON) scripts/build_figures.py

paper:
	bash scripts/build.sh --engine $(ENGINE)

clean:
	-@find $(BUILD_DIR) -type f ! -name "$(JOBNAME).pdf" -delete 2>/dev/null || true
	@echo "Removed intermediate LaTeX files (kept $(PDF) if it exists)."

distclean:
	@rm -rf $(BUILD_DIR)
	@rm -f $(MANUSCRIPT)/figures/*.pdf $(MANUSCRIPT)/figures/*.png
	@rm -f code/figures/*.pdf code/figures/*.png
	@echo "Removed build/ and generated figures."
