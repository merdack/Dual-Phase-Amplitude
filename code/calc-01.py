"""calc-01.py

Numerical illustration of the local (point-dynamics) limit of the
Klein-Gordon-Kuramoto (KKG) model that is presented in the Results section
of the manuscript.

The script integrates a system of N coupled phase-amplitude equations,
computes the Kuramoto order parameter R(t), and produces a three-panel
publication-quality figure that is consumed by the manuscript build
pipeline.

Run as a script:

    python calc-01.py [--prd-width column|full] [--show] [--out PATH] [--format pdf,png]

Defaults:
    - PDF and PNG are written next to this file under ``figures/figure_1.*``
      so they can be copied into ``manuscript/figures`` by the top-level
      build script.
    - The simulation is deterministic (fixed numpy seed) so the figure can
      be regenerated bit-for-bit across machines.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

PrdWidth = Literal["column", "full"]


# ---------------------------------------------------------------------------
# Publication-quality plot defaults
# ---------------------------------------------------------------------------
def prd_figsize_three_panel(width: PrdWidth) -> tuple[float, float]:
    """Approximate PRD widths (inches) for a three-row stacked figure."""
    if width == "column":
        return (3.5, 8.5)
    return (7.0, 10.0)


def configure_matplotlib(prd_width: PrdWidth) -> None:
    """Apply global rcParams; typography scales with ``prd_width``."""
    label = 11 if prd_width == "column" else 14
    tick = 10 if prd_width == "column" else 12
    leg = 9 if prd_width == "column" else 10
    lw = 1.4 if prd_width == "column" else 2.0
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "axes.labelsize": label,
            "font.size": tick,
            "legend.fontsize": leg,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
            "axes.linewidth": 1.2 if prd_width == "column" else 1.5,
            "lines.linewidth": lw,
        }
    )


# ---------------------------------------------------------------------------
# Physics: coupled point-dynamics in the local KKG limit
# ---------------------------------------------------------------------------
def kkg_dynamics(t: float, y: np.ndarray, N: int, K: float, m: np.ndarray) -> np.ndarray:
    """Right-hand side of the simplified KKG point-dynamics ODEs.

    State vector layout:  ``y = [theta_1, ..., theta_N, rho_1, ..., rho_N]``.

    Parameters
    ----------
    t   : current time (unused, kept for ``solve_ivp`` signature).
    y   : flattened phase/amplitude state of length ``2N``.
    N   : number of fields.
    K   : Kuramoto coupling constant.
    m   : per-field intrinsic frequencies (masses), shape ``(N,)``.
    """
    phases = y[:N]
    amps = y[N:]

    sin_diff = np.sin(phases[None, :] - phases[:, None])
    cos_diff = np.cos(phases[None, :] - phases[:, None])

    amp_weighted_sin = (amps[None, :] * sin_diff).sum(axis=1)
    amp_weighted_cos = (amps[None, :] * cos_diff).sum(axis=1)

    d_phases = m + (K / (N * amps)) * amp_weighted_sin
    d_amps = (K / N) * amp_weighted_cos - 0.1 * (amps - 1.0)

    return np.concatenate((d_phases, d_amps))


def simulate(N: int = 5, K: float = 2.5, t_end: float = 20.0, n_points: int = 1000,
             seed: int = 42):
    """Run the simulation and return (t, phases, amplitudes, order parameter)."""
    masses = np.linspace(1.0, 1.5, N)

    rng = np.random.default_rng(seed)
    initial_phases = rng.uniform(-np.pi, np.pi, N)
    initial_amps = rng.uniform(0.8, 1.2, N)
    y0 = np.concatenate((initial_phases, initial_amps))

    t_eval = np.linspace(0.0, t_end, n_points)
    solution = solve_ivp(
        kkg_dynamics,
        (0.0, t_end),
        y0,
        args=(N, K, masses),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )
    if not solution.success:
        raise RuntimeError(f"solve_ivp failed: {solution.message}")

    phases = solution.y[:N, :]
    amps = solution.y[N:, :]
    order = np.abs(np.mean(np.exp(1j * phases), axis=0))
    return t_eval, phases, amps, order


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def make_figure(t, phases, amps, order, *, prd_width: PrdWidth) -> plt.Figure:
    """Build the three-panel figure used in the manuscript Results section."""
    N = phases.shape[0]
    fig, axes = plt.subplots(3, 1, figsize=prd_figsize_three_panel(prd_width), sharex=True)

    for j in range(N):
        axes[0].plot(t, np.sin(phases[j]), alpha=0.8)
    axes[0].set_ylabel(r"$\sin(\theta_j)$")
    axes[0].set_title(r"(a) Phase synchronization of fields ($K > 0$)")
    axes[0].grid(True, linestyle="--", alpha=0.6)

    leg_fs = 8 if prd_width == "column" else 10
    for j in range(N):
        axes[1].plot(t, amps[j], alpha=0.8, label=f"Field {j + 1}")
    axes[1].set_ylabel(r"Amplitude $\rho_j$")
    axes[1].set_title(r"(b) Amplitude evolution and energy exchange")
    axes[1].legend(loc="upper right", ncol=3, fontsize=leg_fs)
    axes[1].grid(True, linestyle="--", alpha=0.6)

    axes[2].plot(t, order, color="black", linewidth=2.5)
    axes[2].set_xlabel(r"Time ($t$)")
    axes[2].set_ylabel(r"Order parameter $R$")
    axes[2].set_title(r"(c) Macroscopic synchronization level")
    axes[2].set_ylim(0, 1.05)
    axes[2].grid(True, linestyle="--", alpha=0.6)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_1",
        help="Output path (without extension). Defaults to code/figures/figure_1.",
    )
    parser.add_argument(
        "--format",
        default="pdf,png",
        help="Comma-separated list of formats to write. Defaults to 'pdf,png'.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively after saving.",
    )
    parser.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD-oriented canvas: 'column' (~3.5 in wide) or 'full' two-column (~7 in).",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"
    configure_matplotlib(prd_width)

    t, phases, amps, order = simulate()
    fig = make_figure(t, phases, amps, order, prd_width=prd_width)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    for fmt in (f.strip() for f in args.format.split(",") if f.strip()):
        out_path = args.out.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"  wrote {out_path}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
