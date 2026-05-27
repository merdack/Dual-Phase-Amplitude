r"""calc-07.py

Three-dimensional comparative visualization of a schematic **probability of presence /
density profile** $\rho(x,t)$ for standard CDM (cusp) versus a Klein--Gordon--Kuramoto
(KKG) motivated **core** from increasing phase stiffness.

* **CDM:** Gaussian whose width shrinks with time, concentrating probability at $x=0$
  (cusp-like sharpening).
* **KKG:** Super-Gaussian $\exp\bigl(-|x/w|^{2m}\bigr)$ with exponent $2m(t)$ tied to a
  proxy order parameter $R(t)\to 1$, producing a progressively flatter central core and
  steeper shoulders (schematic ``synchronization pressure'').

The plot uses two side-by-side ``Axes3D`` surfaces: $x$ (space), $t$ (time), $\rho$
(density).  This is an illustrative cartoon, not a full N-body or field-theory solve.

Run::

    python calc-07.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_7``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

PrdWidth = Literal["column", "full"]


def prd_figsize_dual_3d(width: PrdWidth) -> tuple[float, float]:
    """PRD-oriented figure size (inches) for two side-by-side 3D axes."""
    if width == "column":
        return (7.4, 4.8)
    return (14.0, 6.2)


def configure_matplotlib(prd_width: PrdWidth) -> None:
    label = 11 if prd_width == "column" else 13
    tick = 9 if prd_width == "column" else 11
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "axes.labelsize": label,
            "font.size": tick,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
            "axes.linewidth": 1.0,
        }
    )


def build_grids(
    x_max: float,
    nx: int,
    t_max: float,
    nt: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.linspace(-float(x_max), float(x_max), int(nx))
    t = np.linspace(0.0, float(t_max), int(nt))
    X, T = np.meshgrid(x, t)
    return x, t, X, T


def density_cdm(X: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Narrowing Gaussian: probability piles up at the origin (cusp-like)."""
    sigma = 0.48 * np.exp(-0.62 * T) + 0.022
    norm = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    return norm * np.exp(-0.5 * (X / sigma) ** 2)


def order_parameter_proxy(T: np.ndarray, tau: float = 1.05) -> np.ndarray:
    """Smooth proxy $R(t)\\to 1$ (Kuramoto-like synchronization level)."""
    return 1.0 - np.exp(-T / tau)


def density_kkg(X: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Super-Gaussian core: exponent grows with $R(t)$ (stiffer effective fluid)."""
    R = order_parameter_proxy(T)
    m = 1.0 + 4.2 * R**2
    w = 0.52 + 0.22 * R
    exponent = 2.0 * m
    return np.exp(-(np.abs(X) / np.maximum(w, 1e-6)) ** exponent)


def make_figure(
    X: np.ndarray,
    T: np.ndarray,
    z_cdm: np.ndarray,
    z_kkg: np.ndarray,
    *,
    prd_width: PrdWidth,
    elev: float,
    azim: float,
    stride: int,
) -> plt.Figure:
    fig = plt.figure(figsize=prd_figsize_dual_3d(prd_width))
    ax_cdm = fig.add_subplot(1, 2, 1, projection="3d")
    ax_kkg = fig.add_subplot(1, 2, 2, projection="3d")

    rs = max(1, int(stride))
    surf_cdm = ax_cdm.plot_surface(
        X,
        T,
        z_cdm,
        cmap="viridis",
        linewidth=0,
        antialiased=True,
        rstride=rs,
        cstride=rs,
        edgecolor="none",
    )
    surf_kkg = ax_kkg.plot_surface(
        X,
        T,
        z_kkg,
        cmap="magma",
        linewidth=0,
        antialiased=True,
        rstride=rs,
        cstride=rs,
        edgecolor="none",
    )

    for ax, title in (
        (
            ax_cdm,
            "CDM (schematic): cusp / central concentration",
        ),
        (
            ax_kkg,
            r"PLKG (schematic): core via $R(t)\!\to\!1$ (super-Gaussian stiffening)",
        ),
    ):
        ax.set_title(title, pad=12)
        ax.set_xlabel(r"Space $x$")
        ax.set_ylabel(r"Time $t$")
        ax.set_zlabel(r"Density $\rho(x,t)$")
        ax.view_init(elev=elev, azim=azim)

    fig.colorbar(surf_cdm, ax=ax_cdm, shrink=0.55, pad=0.12, label=r"$\rho$ (CDM)")
    fig.colorbar(surf_kkg, ax=ax_kkg, shrink=0.55, pad=0.12, label=r"$\rho$ (PLKG)")

    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.06, top=0.9, wspace=0.28)
    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Schematic 3D CDM (cusp) vs KKG (core) density surfaces.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_7",
        help="Output path without extension (default: code/figures/figure_7).",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD canvas width hint for figure size (column ~7.4 in, full ~14 in).",
    )
    p.add_argument("--x-max", type=float, default=2.5, help="Half-width of spatial domain.")
    p.add_argument("--t-max", type=float, default=2.6, help="Maximum time.")
    p.add_argument("--nx", type=int, default=90, help="Spatial grid points.")
    p.add_argument("--nt", type=int, default=70, help="Temporal grid points.")
    p.add_argument(
        "--stride",
        type=int,
        default=2,
        help="plot_surface rstride/cstride (>=1; larger is faster/coarser).",
    )
    p.add_argument("--elev", type=float, default=28.0, help="3D view elevation (deg).")
    p.add_argument("--azim", type=float, default=-58.0, help="3D view azimuth (deg).")
    p.add_argument("--show", action="store_true", help="Display interactively.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    if args.nx < 4 or args.nt < 4:
        raise SystemExit("error: require nx, nt >= 4")
    if args.stride < 1:
        raise SystemExit("error: --stride must be >= 1")

    configure_matplotlib(prd_width)
    _x, _t, X, T = build_grids(args.x_max, args.nx, args.t_max, args.nt)
    z_cdm = density_cdm(X, T)
    z_kkg = density_kkg(X, T)

    fig = make_figure(
        X,
        T,
        z_cdm,
        z_kkg,
        prd_width=prd_width,
        elev=args.elev,
        azim=args.azim,
        stride=args.stride,
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
