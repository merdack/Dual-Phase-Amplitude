"""calc-04.py

Pedagogical cosmology plot for the synchronization efficiency scale
$K_{\\mathrm{cosmo}}\\sim 3H^2N$ (manuscript App.~B; Eq.~sync_threshold).

Panel (a): $H(z)/H_0$ for a spatially flat FRW background.
Panel (b): $K/\\bigl(3H(z)^2N/H_0^2\\bigr)$ with $H(z)=H_0 E(z)$ and
$E(z)=\\sqrt{\\Omega_m(1+z)^3+\\Omega_\\Lambda}$, treating each curve as a constant $K/H_0^2$.
The horizontal dotted line marks $K=3H(z)^2N$ in these schematic units.

Run::

    python calc-04.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_4``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

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


def e_of_z(z: np.ndarray, omega_m: float, omega_l: float) -> np.ndarray:
    """Flat FRW: $E(z)=H(z)/H_0$ with matter + cosmological constant."""
    return np.sqrt(np.maximum(omega_m * (1.0 + z) ** 3 + omega_l, 1e-30))


def k_over_three_h2n(
    z: np.ndarray,
    k_units: float,
    *,
    n: int,
    omega_m: float,
    omega_l: float,
) -> np.ndarray:
    """Return $K / (3 N H(z)^2)$ in schematic units with $[K]=H_0^2$.

    We write $H(z)=H_0 E(z)$ and treat the numerator as a constant $K$ measured
    in units of $H_0^2$, so $K/(3 N H(z)^2)= (K/H_0^2)/(3 N E(z)^2)$ is dimensionless.
    """
    e = e_of_z(z, omega_m, omega_l)
    return k_units / (3.0 * float(n) * e * e)


def make_figure(
    z: np.ndarray,
    *,
    ks: tuple[float, ...],
    n: int,
    omega_m: float,
    omega_l: float,
    prd_width: PrdWidth,
) -> plt.Figure:
    fig, (ax0, ax1) = plt.subplots(
        2,
        1,
        figsize=prd_figsize_two_panel(prd_width),
        sharex=True,
        layout="constrained",
    )

    e = e_of_z(z, omega_m, omega_l)
    ax0.plot(z, e, color="black")
    ax0.set_ylabel(r"$H(z)/H_0$")
    ax0.set_title(r"(a) Background expansion (flat FRW)")
    ax0.grid(True, linestyle="--", alpha=0.55)

    cmap = plt.get_cmap("tab10")
    for i, k_u in enumerate(ks):
        ratio = k_over_three_h2n(z, k_u, n=n, omega_m=omega_m, omega_l=omega_l)
        ax1.plot(z, ratio, color=cmap(i % 10), label=rf"$K/H_0^2={k_u:g}$")

    ax1.axhline(1.0, color="0.35", linestyle=":", linewidth=1.2, label=r"$K=3H(z)^2N$")
    ax1.set_xlabel(r"Redshift $z$")
    ax1.set_ylabel(r"$\frac{K}{3H(z)^2N/H_0^2}$")
    ax1.set_title(rf"(b) Cosmological efficiency proxy ($N={n}$)")
    ax1.grid(True, linestyle="--", alpha=0.55)
    ax1.legend(loc="upper right", framealpha=0.95)
    ax1.set_xlim(float(z[0]), float(z[-1]))

    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_4",
        help="Output path without extension (default: code/figures/figure_4).",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD canvas width: 'column' (~3.5 in) or 'full' (~7.0 in).",
    )
    p.add_argument("--omega-m", type=float, default=0.315, help="Matter density parameter (default 0.315).")
    p.add_argument("--omega-l", type=float, default=0.685, help="Dark energy density parameter (default 0.685).")
    p.add_argument("--z-max", type=float, default=45.0, help="Maximum redshift on the grid.")
    p.add_argument("--n-z", type=int, default=400, help="Number of redshift samples.")
    p.add_argument("--N", type=int, default=5, dest="N_fields", help="Number of fields $N$ (default 5).")
    p.add_argument(
        "--K",
        default="10,35,120",
        help="Comma-separated list of $K/H_0^2$ values for panel (b) (default: 10,35,120).",
    )
    p.add_argument("--show", action="store_true", help="Show interactively.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    if abs(args.omega_m + args.omega_l - 1.0) > 1e-3:
        raise SystemExit("error: require Omega_m + Omega_Lambda ~= 1 (flatness) for this schematic")
    if args.n_z < 2 or args.z_max <= 0.0:
        raise SystemExit("error: invalid z grid")
    if args.N_fields < 1:
        raise SystemExit("error: N must be positive")

    ks = tuple(float(x.strip()) for x in args.K.split(",") if x.strip())
    if not ks:
        raise SystemExit("error: --K must list at least one value")

    z = np.linspace(0.0, args.z_max, args.n_z)
    configure_matplotlib(prd_width)
    fig = make_figure(
        z,
        ks=ks,
        n=args.N_fields,
        omega_m=args.omega_m,
        omega_l=args.omega_l,
        prd_width=prd_width,
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
