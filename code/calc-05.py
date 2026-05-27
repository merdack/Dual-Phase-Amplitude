"""calc-05.py

Schematic linear-stability figure (roadmap idea 4): effective mass squared
$m_{\\mathrm{eff}}^2=m^2-K-\\dot{\\bar{\\theta}}^2$ from App.~B and the Jeans-type
instability boundary $m_{\\mathrm{eff}}^2=0$.

Panel (a): $m_{\\mathrm{eff}}^2$ versus $K$ for several fixed values of $\\dot{\\bar{\\theta}}^2$.
Panel (b): $(K,\\dot{\\theta}^2)$ plane with a colormap of $m_{\\mathrm{eff}}^2$ and a
black contour at zero (App.~B: ``effective mass'' and Jeans-type instability).

Run::

    python calc-05.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_5``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm

PrdWidth = Literal["column", "full"]


def prd_figsize_two_panel(width: PrdWidth) -> tuple[float, float]:
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
            "legend.fontsize": tick - 1,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
            "axes.linewidth": 1.2,
            "lines.linewidth": lw,
        }
    )


def make_figure(
    k: np.ndarray,
    *,
    m: float,
    theta_dot_sqs: tuple[float, ...],
    k2: np.ndarray,
    td2: np.ndarray,
    prd_width: PrdWidth,
) -> plt.Figure:
    fig, (ax0, ax1) = plt.subplots(
        2,
        1,
        figsize=prd_figsize_two_panel(prd_width),
        layout="constrained",
    )

    cmap_lines = plt.get_cmap("tab10")
    for i, td2_i in enumerate(theta_dot_sqs):
        meff2 = m * m - k - td2_i
        ax0.plot(k, meff2, color=cmap_lines(i % 10), label=rf"$\dot{{\bar{{\theta}}}}^2={td2_i:g}$")
        k_crit = m * m - td2_i
        if 0.0 < k_crit < float(k[-1]):
            ax0.axvline(k_crit, color=cmap_lines(i % 10), linestyle=":", linewidth=1.0, alpha=0.85)

    ax0.axhline(0.0, color="0.35", linestyle="--", linewidth=1.1)
    ax0.set_xlabel(r"Coupling $K$")
    ax0.set_ylabel(r"$m_{\mathrm{eff}}^2$")
    ax0.set_title(r"(a) Effective mass squared vs $K$ (fixed $\dot{\bar{\theta}}^2$)")
    ax0.grid(True, linestyle="--", alpha=0.55)
    ax0.legend(loc="best", framealpha=0.95)

    kk, tt = np.meshgrid(k2, td2)
    meff2_map = m * m - kk - tt
    vmin = float(np.min(meff2_map))
    vmax = float(np.max(meff2_map))
    if vmin < 0.0 < vmax:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=max(vmax, 1e-12))
    cf = ax1.contourf(kk, tt, meff2_map, levels=40, cmap="RdBu_r", norm=norm)
    ax1.contour(kk, tt, meff2_map, levels=[0.0], colors="k", linewidths=1.3)
    ax1.set_xlabel(r"Coupling $K$")
    ax1.set_ylabel(r"$\dot{\bar{\theta}}^2$")
    ax1.set_title(r"(b) $(K,\dot{\bar{\theta}}^2)$ map of $m_{\mathrm{eff}}^2$ ($m=1$ units)")
    cbar = fig.colorbar(cf, ax=ax1, shrink=0.85, pad=0.02)
    cbar.set_label(r"$m_{\mathrm{eff}}^2$")

    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_5",
        help="Output path without extension (default: code/figures/figure_5).",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD canvas width: 'column' (~3.5 in) or 'full' (~7.0 in).",
    )
    p.add_argument("--m", type=float, default=1.0, help="Background mass scale m (default 1).")
    p.add_argument(
        "--theta-dotsq",
        default="0,0.1,0.25",
        help="Comma-separated dot(theta)^2 values for panel (a).",
    )
    p.add_argument("--k-min", type=float, default=0.0, help="Minimum K in panel (a).")
    p.add_argument("--k-max", type=float, default=1.35, help="Maximum K in panel (a).")
    p.add_argument("--n-k", type=int, default=300, help="Number of K samples in panel (a).")
    p.add_argument(
        "--map-k-max",
        type=float,
        default=None,
        help="K extent for panel (b); defaults to --k-max.",
    )
    p.add_argument(
        "--map-td2-max",
        type=float,
        default=None,
        help=r"Maximum $\dot{\theta}^2$ on panel (b); defaults to m^2.",
    )
    p.add_argument("--n-map", type=int, default=180, help="Grid size per axis for panel (b).")
    p.add_argument("--show", action="store_true", help="Show interactively.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    if args.m <= 0.0:
        raise SystemExit("error: require m > 0")
    if args.k_max <= args.k_min:
        raise SystemExit("error: require --k-max > --k-min")
    if args.n_k < 2 or args.n_map < 2:
        raise SystemExit("error: invalid grid resolution")

    theta_dot_sqs = tuple(float(x.strip()) for x in args.theta_dotsq.split(",") if x.strip())
    if not theta_dot_sqs:
        raise SystemExit("error: --theta-dotsq must list at least one value")

    k = np.linspace(args.k_min, args.k_max, args.n_k)
    k2_max = float(args.map_k_max) if args.map_k_max is not None else float(args.k_max)
    td2_max = float(args.map_td2_max) if args.map_td2_max is not None else float(args.m) ** 2
    k2 = np.linspace(0.0, k2_max, args.n_map)
    td2 = np.linspace(0.0, td2_max, args.n_map)

    configure_matplotlib(prd_width)
    fig = make_figure(k, m=args.m, theta_dot_sqs=theta_dot_sqs, k2=k2, td2=td2, prd_width=prd_width)

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
