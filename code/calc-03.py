"""calc-03.py

Coupling-resolved synchronization summary for the local KKG point model
(roadmap idea 2): asymptotic order parameter and first crossing time versus K,
using the same initial draw and homogeneous ODE truncation as Figs.~1--2.

Run::

    python calc-03.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_3``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

PrdWidth = Literal["column", "full"]


def prd_figsize_two_panel(width: PrdWidth) -> tuple[float, float]:
    """PRD-oriented width with room for two stacked panels (inches)."""
    w_in = 3.5 if width == "column" else 7.0
    h_in = 5.2 if width == "column" else 5.6
    return (w_in, h_in)


def configure_matplotlib(prd_width: PrdWidth) -> None:
    label = 11 if prd_width == "column" else 14
    tick = 10 if prd_width == "column" else 12
    lw = 1.4 if prd_width == "column" else 2.0
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "axes.labelsize": label,
            "font.size": tick,
            "legend.fontsize": tick,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
            "axes.linewidth": 1.2,
            "lines.linewidth": lw,
        }
    )


def kkg_dynamics(t: float, y: np.ndarray, N: int, K: float, m: np.ndarray) -> np.ndarray:
    """RHS of the simplified KKG point-dynamics ODEs (matches ``calc-01.py``)."""
    phases = y[:N]
    amps = y[N:]

    sin_diff = np.sin(phases[None, :] - phases[:, None])
    cos_diff = np.cos(phases[None, :] - phases[:, None])

    amp_weighted_sin = (amps[None, :] * sin_diff).sum(axis=1)
    amp_weighted_cos = (amps[None, :] * cos_diff).sum(axis=1)

    d_phases = m + (K / (N * amps)) * amp_weighted_sin
    d_amps = (K / N) * amp_weighted_cos - 0.1 * (amps - 1.0)

    return np.concatenate((d_phases, d_amps))


def simulate(
    N: int = 5,
    K: float = 2.5,
    t_end: float = 20.0,
    n_points: int = 1000,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(t, order_parameter)`` for one coupling ``K``."""
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
        raise RuntimeError(f"solve_ivp failed for K={K}: {solution.message}")

    phases = solution.y[:N, :]
    order = np.abs(np.mean(np.exp(1j * phases), axis=0))
    return t_eval, order


def first_crossing_time(t: np.ndarray, order: np.ndarray, threshold: float) -> float:
    """First ``t`` where ``order >= threshold``; ``nan`` if never crossed."""
    above = np.nonzero(order >= threshold)[0]
    if above.size == 0:
        return float("nan")
    j = int(above[0])
    if j == 0:
        return float(t[0])
    # Linear interpolation between samples for a smoother crossing estimate
    j0, j1 = j - 1, j
    r0, r1 = float(order[j0]), float(order[j1])
    if r1 <= r0:
        return float(t[j])
    frac = (threshold - r0) / (r1 - r0)
    return float(t[j0] + frac * (t[j1] - t[j0]))


def scan_k(
    k_values: np.ndarray,
    *,
    t_end: float,
    n_points: int,
    seed: int,
    r_threshold: float,
    N: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(K, R_final, t_cross)`` with one run per ``K`` (same seed each time)."""
    r_final = np.empty_like(k_values, dtype=float)
    t_cross = np.empty_like(k_values, dtype=float)
    for i, K in enumerate(k_values):
        t, order = simulate(N=N, K=float(K), t_end=t_end, n_points=n_points, seed=seed)
        r_final[i] = float(order[-1])
        t_cross[i] = first_crossing_time(t, order, r_threshold)
    return k_values, r_final, t_cross


def make_figure(
    k: np.ndarray,
    r_final: np.ndarray,
    t_cross: np.ndarray,
    *,
    prd_width: PrdWidth,
    r_threshold: float,
    t_end: float,
) -> plt.Figure:
    fig, (ax0, ax1) = plt.subplots(
        2,
        1,
        figsize=prd_figsize_two_panel(prd_width),
        sharex=True,
        layout="constrained",
    )

    ax0.plot(k, r_final, color="C0")
    ax0.set_ylabel(r"$R(t_{\mathrm{f}})$")
    ax0.set_ylim(0.0, 1.05)
    ax0.grid(True, linestyle="--", alpha=0.55)
    ax0.set_title(rf"(a) Terminal order parameter at $t_{{\mathrm{{f}}}}={t_end:g}$")

    mask = np.isfinite(t_cross)
    ax1.plot(k[mask], t_cross[mask], color="C1")
    ax1.set_xlabel(r"Coupling $K$")
    ax1.set_ylabel(rf"$t_{{{r_threshold:g}}}$")
    ax1.grid(True, linestyle="--", alpha=0.55)
    ax1.set_title(rf"(b) First time $R(t)\geq {r_threshold:g}$ (missing if not reached)")

    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_3",
        help="Output path without extension (default: code/figures/figure_3).",
    )
    p.add_argument(
        "--format",
        default="pdf,png",
        help="Comma-separated formats (default: pdf,png).",
    )
    p.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD canvas width: 'column' (~3.5 in) or 'full' (~7.0 in).",
    )
    p.add_argument("--k-min", type=float, default=0.05, help="Minimum K in scan.")
    p.add_argument("--k-max", type=float, default=12.0, help="Maximum K in scan.")
    p.add_argument("--n-k", type=int, default=56, help="Number of K samples (uniform grid).")
    p.add_argument(
        "--t-end",
        type=float,
        default=45.0,
        help="Fixed integration horizon for every K (default: 45).",
    )
    p.add_argument(
        "--n-points",
        type=int,
        default=1200,
        help="Number of uniform time samples (default: 1200).",
    )
    p.add_argument(
        "--r-threshold",
        type=float,
        default=0.9,
        help="Threshold for first-crossing time panel (default: 0.9).",
    )
    p.add_argument("--seed", type=int, default=42, help="RNG seed for initial data (default: 42).")
    p.add_argument("--show", action="store_true", help="Show the figure interactively.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    if args.n_k < 2:
        raise SystemExit("error: --n-k must be at least 2")
    if args.k_max <= args.k_min:
        raise SystemExit("error: require --k-max > --k-min")

    k_values = np.linspace(args.k_min, args.k_max, args.n_k)
    configure_matplotlib(prd_width)
    k_arr, r_final, t_cross = scan_k(
        k_values,
        t_end=args.t_end,
        n_points=args.n_points,
        seed=args.seed,
        r_threshold=args.r_threshold,
    )
    fig = make_figure(
        k_arr,
        r_final,
        t_cross,
        prd_width=prd_width,
        r_threshold=args.r_threshold,
        t_end=args.t_end,
    )

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
