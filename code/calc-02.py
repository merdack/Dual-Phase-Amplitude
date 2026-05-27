"""calc-02.py

Order parameter R(t) versus time for several Kuramoto couplings K, using the
same local KKG point-dynamics model and identical initial data for each curve
(roadmap idea 1: coupling-strength scan of synchronization).

Run::

    python calc-02.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default output basename is ``figures/figure_2``.
The ``--prd-width`` flag sets approximate Physical Review column geometry in inches
before the manuscript scales ``\\includegraphics``.
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


def prd_figsize(width: PrdWidth, *, aspect: float = 0.72) -> tuple[float, float]:
    """Approximate PRD single-column vs two-column widths (inches)."""
    w_in = 3.5 if width == "column" else 7.0
    return (w_in, max(w_in * aspect, 2.4))


def configure_matplotlib(prd_width: PrdWidth) -> None:
    """Global rcParams; slightly smaller typography in single-column mode."""
    label = 11 if prd_width == "column" else 14
    tick = 10 if prd_width == "column" else 12
    leg = 9 if prd_width == "column" else 11
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(t, phases, amplitudes, order)`` for one coupling ``K``."""
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


def simulate_order_curves(
    ks: tuple[float, ...],
    *,
    N: int = 5,
    t_end: float = 20.0,
    n_points: int = 1000,
    seed: int = 42,
) -> tuple[np.ndarray, dict[float, np.ndarray]]:
    """Run one integration per K with identical initial data (same RNG seed)."""
    t_ref: np.ndarray | None = None
    orders: dict[float, np.ndarray] = {}
    for K in ks:
        t, _ph, _am, order = simulate(N=N, K=K, t_end=t_end, n_points=n_points, seed=seed)
        if t_ref is None:
            t_ref = t
        elif not np.allclose(t_ref, t):
            raise RuntimeError("internal: t grids differ between runs")
        orders[K] = order
    assert t_ref is not None
    return t_ref, orders


def make_figure(
    t: np.ndarray,
    orders: dict[float, np.ndarray],
    *,
    prd_width: PrdWidth,
) -> plt.Figure:
    """Single-panel R(t) for several K."""
    figsize = prd_figsize(prd_width, aspect=0.62)
    fig, ax = plt.subplots(figsize=figsize)

    cmap = plt.get_cmap("tab10")
    for i, (K, R) in enumerate(sorted(orders.items())):
        ax.plot(t, R, color=cmap(i % 10), label=rf"$K={K:g}$")

    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Order parameter $R(t)$")
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(float(t[0]), float(t[-1]))
    ax.grid(True, linestyle="--", alpha=0.55)
    ax.legend(loc="lower right", framealpha=0.92)
    fig.tight_layout()
    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_2",
        help="Output path without extension (default: code/figures/figure_2).",
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
        help="PRD-oriented canvas width in inches: 'column' (~3.5 in) or "
        "'full' two-column (~7.0 in). Height scales with width.",
    )
    p.add_argument(
        "--K",
        default="0.2,0.8,2.5,8.0",
        help="Comma-separated list of coupling constants (default: 0.2,0.8,2.5,8.0).",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive window after saving.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    ks = tuple(float(x.strip()) for x in args.K.split(",") if x.strip())
    if not ks:
        raise SystemExit("error: --K must contain at least one value")

    configure_matplotlib(prd_width)
    t, orders = simulate_order_curves(ks)
    fig = make_figure(t, orders, prd_width=prd_width)

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
