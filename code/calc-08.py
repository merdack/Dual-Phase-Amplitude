r"""calc-08.py

Toy microscopic simulation comparing **probability density** $|\Phi(x,t)|^2$ for a
multi-mode field $\Phi=\sum_j \rho_j(x)\,e^{i\theta_j(t)}$ under:

* **FDM-like sector ($K=0$):** independent phase drift plus common white-noise perturbations
  (no Kuramoto coupling).
* **KKG sector ($K>0$):** the same noise draws plus a Kuramoto synchronizing torque
  $(K/N)\sum_k\sin(\theta_k-\theta_j)$.

The script builds fixed Gaussian envelopes $\rho_j(x)$ on $x\in[-10,10]$, advances
phases with Euler--Maruyama time steps, forms $\Phi$ on a spatial grid, and plots
side-by-side $|\Phi|^2$ heat maps together with $\max_x|\Phi|^2$ versus time and the
temporal variance of that maximum (post burn-in).

This is a **schematic reduced model** for visualization, not a full quantum FDM or
covariant KKG field solve.

Run::

    python calc-08.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_8``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

PrdWidth = Literal["column", "full"]


def prd_figsize(width: PrdWidth) -> tuple[float, float]:
    if width == "column":
        return (7.4, 9.0)
    return (14.0, 9.5)


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
            "legend.fontsize": tick,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
        }
    )


def gaussian_envelopes(
    x: np.ndarray,
    centers: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Return ``G[x, j]`` with shape ``(nx, N)``, Gaussian packets centered at ``centers[j]``."""
    dx = x[:, None] - centers[None, :]
    return np.exp(-0.5 * (dx / sigma) ** 2, dtype=np.float64)


def simulate_phi_squared(
    *,
    K: float,
    noise: np.ndarray,
    omega: np.ndarray,
    G: np.ndarray,
    dt: float,
    eta: float,
    theta0: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(rho_hist, M_max)`` with ``rho_hist`` shape ``(n_steps+1, nx)``."""
    n_steps, _ = noise.shape
    nx = G.shape[0]
    theta = theta0.astype(np.float64, copy=True)
    rho_hist = np.empty((n_steps + 1, nx), dtype=np.float64)
    M_max = np.empty(n_steps + 1, dtype=np.float64)
    sqrt_dt = np.sqrt(dt)

    for idx in range(n_steps + 1):
        phi = G @ np.exp(1j * theta)
        rho = np.abs(phi) ** 2
        rho_hist[idx] = rho
        M_max[idx] = float(np.max(rho))
        if idx == n_steps:
            break
        xi = noise[idx]
        kuramoto = (K / G.shape[1]) * np.sin(theta[:, None] - theta[None, :]).sum(axis=0)
        theta = theta + omega * dt + eta * sqrt_dt * xi + kuramoto * dt

    return rho_hist, M_max


def temporal_variance_max_density(M: np.ndarray, burn_in: int) -> float:
    """Variance of ``M[t]`` after discarding the first ``burn_in`` samples."""
    if burn_in >= len(M):
        return float("nan")
    return float(np.var(M[burn_in:], ddof=0))


def make_figure(
    x: np.ndarray,
    t: np.ndarray,
    rho_fdm: np.ndarray,
    rho_kkg: np.ndarray,
    M_fdm: np.ndarray,
    M_kkg: np.ndarray,
    var_fdm: float,
    var_kkg: float,
    *,
    prd_width: PrdWidth,
    burn_in: int,
) -> plt.Figure:
    fig = plt.figure(figsize=prd_figsize(prd_width))
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[2.3, 1.0], hspace=0.35, wspace=0.28)

    extent = (float(x[0]), float(x[-1]), float(t[0]), float(t[-1]))
    vmax = float(np.percentile(np.concatenate((rho_fdm.ravel(), rho_kkg.ravel())), 99.5))
    vmax = max(vmax, 1e-12)
    vmin = 0.0

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    im0 = ax_a.imshow(
        rho_fdm,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    im1 = ax_b.imshow(
        rho_kkg,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    ax_a.set_title(r"(A) FDM-like ($K=0$): incoherent phases + noise")
    ax_b.set_title(r"(B) PLKG-like ($K>0$): phase-locked coupling + same noise")
    ax_a.set_xlabel(r"Space $x$")
    ax_b.set_xlabel(r"Space $x$")
    ax_a.set_ylabel(r"Time $t$")
    ax_b.set_ylabel(r"Time $t$")
    fig.colorbar(im0, ax=[ax_a, ax_b], shrink=0.72, pad=0.03, label=r"$|\Phi|^2$")

    ax_c = fig.add_subplot(gs[1, :])
    ax_c.plot(
        t,
        M_fdm,
        color="C0",
        label=rf"FDM: temporal Var of $\max|\Phi|^2$ (steps $\geq$ {burn_in}) = {var_fdm:.4g}",
    )
    ax_c.plot(
        t,
        M_kkg,
        color="C3",
        label=rf"PLKG: temporal Var of $\max|\Phi|^2$ (steps $\geq$ {burn_in}) = {var_kkg:.4g}",
    )
    ax_c.set_xlabel(r"Time $t$")
    ax_c.set_ylabel(r"$M(t)=\max_x |\Phi|^2$")
    ax_c.set_title("(C) Maximum density vs time (lower temporal Var(M) implies stabler peak)")
    ax_c.grid(True, linestyle="--", alpha=0.55)
    ax_c.legend(loc="upper right", framealpha=0.95)

    fig.subplots_adjust(top=0.93, bottom=0.07, left=0.06, right=0.97)
    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="FDM vs KKG toy simulation of |Phi|^2 on a line with shared white noise.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_8",
        help="Output path without extension (default: code/figures/figure_8).",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD-oriented figure width (column ~7.4 in, full ~14 in).",
    )
    p.add_argument("--N", type=int, default=40, help="Number of modes / packets.")
    p.add_argument("--m", type=float, default=1.0, help="Base intrinsic frequency scale.")
    p.add_argument("--K-kkg", type=float, default=5.0, dest="K_kkg", help="Kuramoto coupling for KKG panel.")
    p.add_argument("--eta", type=float, default=0.3, help="Noise strength on phases.")
    p.add_argument("--dt", type=float, default=0.05, help="Time step.")
    p.add_argument("--n-steps", type=int, default=200, help="Number of Euler--Maruyama steps.")
    p.add_argument("--x-min", type=float, default=-10.0)
    p.add_argument("--x-max", type=float, default=10.0)
    p.add_argument("--nx", type=int, default=480, help="Spatial grid resolution.")
    p.add_argument("--sigma", type=float, default=0.48, help="Gaussian packet width.")
    p.add_argument("--burn-in", type=int, default=30, help="Steps dropped for Var(M) estimate.")
    p.add_argument("--seed", type=int, default=2026, help="RNG seed (reproducible builds).")
    p.add_argument("--show", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    if args.N < 2 or args.nx < 16 or args.n_steps < 1:
        raise SystemExit("error: invalid N, nx, or n-steps")
    if args.x_max <= args.x_min:
        raise SystemExit("error: require x_max > x_min")

    rng = np.random.default_rng(args.seed)
    x = np.linspace(args.x_min, args.x_max, args.nx, dtype=np.float64)
    span = args.x_max - args.x_min
    margin = 0.05 * span
    centers = np.linspace(args.x_min + margin, args.x_max - margin, args.N, dtype=np.float64)
    G = gaussian_envelopes(x, centers, args.sigma)

    omega = args.m * (1.0 + 0.03 * rng.standard_normal(args.N))
    theta0 = rng.uniform(-np.pi, np.pi, size=args.N).astype(np.float64)
    noise = rng.standard_normal((args.n_steps, args.N)).astype(np.float64)

    rho_fdm, M_fdm = simulate_phi_squared(
        K=0.0,
        noise=noise,
        omega=omega,
        G=G,
        dt=args.dt,
        eta=args.eta,
        theta0=theta0,
    )
    rho_kkg, M_kkg = simulate_phi_squared(
        K=args.K_kkg,
        noise=noise,
        omega=omega,
        G=G,
        dt=args.dt,
        eta=args.eta,
        theta0=theta0.copy(),
    )

    t = np.linspace(0.0, args.dt * args.n_steps, args.n_steps + 1)
    var_fdm = temporal_variance_max_density(M_fdm, args.burn_in)
    var_kkg = temporal_variance_max_density(M_kkg, args.burn_in)

    configure_matplotlib(prd_width)
    fig = make_figure(
        x,
        t,
        rho_fdm,
        rho_kkg,
        M_fdm,
        M_kkg,
        var_fdm,
        var_kkg,
        prd_width=prd_width,
        burn_in=args.burn_in,
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

    print(f"  Var(M) FDM (post burn-in) = {var_fdm:.6g}")
    print(f"  Var(M) KKG (post burn-in) = {var_kkg:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
