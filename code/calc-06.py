r"""calc-06.py

Stability-focused summary for the local PLKG point model (Fig.~6): mean amplitude
versus $K$ in the classically stable window $K<m^2$, the linear-theory effective
mass squared $m_{\mathrm{eff}}^2$ from Appendix~B, and a cubic-regularized
amplitude trajectory at strong coupling ($K=8$) compared to the $\lambda=0$
toy truncation.

Run::

    python calc-06.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_6``.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

PrdWidth = Literal["column", "full"]


def prd_figsize_three_panel(width: PrdWidth) -> tuple[float, float]:
    w_in = 3.5 if width == "column" else 7.0
    h_in = 6.15 if width == "column" else 6.65
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


def kkg_dynamics(
    t: float,
    y: np.ndarray,
    N: int,
    K: float,
    m: np.ndarray,
    lam: float,
) -> np.ndarray:
    """RHS of the KKG point ODEs with optional quartic amplitude self-interaction."""
    phases = y[:N]
    amps = y[N:]

    sin_diff = np.sin(phases[None, :] - phases[:, None])
    cos_diff = np.cos(phases[None, :] - phases[:, None])

    amp_weighted_sin = (amps[None, :] * sin_diff).sum(axis=1)
    amp_weighted_cos = (amps[None, :] * cos_diff).sum(axis=1)

    d_phases = m + (K / (N * amps)) * amp_weighted_sin
    d_amps = (K / N) * amp_weighted_cos - 0.1 * (amps - 1.0) - lam * (amps - 1.0) ** 3
    return np.concatenate((d_phases, d_amps))


def simulate(
    N: int,
    K: float,
    masses: np.ndarray,
    *,
    t_end: float,
    n_points: int,
    seed: int,
    lam: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Return ``(t, amps, ok)`` with ``amps`` shape ``(N, len(t))``.

    For panel~(c), trajectories are truncated at the first sample where the mean
    amplitude or any $\rho_j$ leaves a moderate band (numerical blow-up guard).
    """
    rng = np.random.default_rng(seed)
    initial_phases = rng.uniform(-np.pi, np.pi, N)
    initial_amps = rng.uniform(0.8, 1.2, N)
    y0 = np.concatenate((initial_phases, initial_amps))

    t_eval = np.linspace(0.0, t_end, n_points)
    solution = solve_ivp(
        kkg_dynamics,
        (0.0, t_end),
        y0,
        args=(N, K, masses, lam),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )
    if solution.t.size == 0:
        raise RuntimeError(f"solve_ivp returned no samples for K={K}, lambda={lam}")

    amps = solution.y[N:, :]
    t = np.asarray(solution.t, dtype=np.float64)

    mean_r = np.mean(amps, axis=0)
    blow = (
        ~np.isfinite(mean_r)
        | np.any(~np.isfinite(amps), axis=0)
        | np.any(amps < 5e-3, axis=0)
        | np.any(amps > 80.0, axis=0)
        | (mean_r > 80.0)
    )
    if np.any(blow):
        first = int(np.argmax(blow))
        if first <= 2:
            t = t[:1]
            amps = amps[:, :1]
        else:
            t = t[:first]
            amps = amps[:, :first]

    ok = bool(solution.success) and np.all(np.isfinite(amps)) and t.size >= 4
    return t, amps, ok


def mean_rho_tail(amps: np.ndarray, burn_fraction: float) -> float:
    """Mean $\rho$ over oscillators and over the post-burn tail."""
    if not (0.0 <= burn_fraction < 1.0):
        raise ValueError("burn_fraction must be in [0, 1)")
    n_t = amps.shape[1]
    tb = int(np.floor(n_t * burn_fraction))
    tb = min(max(tb, 0), n_t - 1)
    tail = np.asarray(amps[:, tb:], dtype=np.float64)
    if not np.all(np.isfinite(tail)):
        raise ValueError("non-finite amplitudes in tail window")
    return float(np.mean(tail))


def scan_mean_rho_stable(
    k_values: np.ndarray,
    *,
    N: int,
    masses: np.ndarray,
    t_end: float,
    n_points: int,
    seed: int,
    burn_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Mean tail $\langle\rho\rangle$ for each ``K`` (identical masses and seed)."""
    mean_rho = np.empty_like(k_values, dtype=float)
    for i, K in enumerate(k_values):
        _t, amps, _ok = simulate(
            N,
            float(K),
            masses,
            t_end=t_end,
            n_points=n_points,
            seed=seed,
            lam=0.0,
        )
        mean_rho[i] = mean_rho_tail(amps, burn_fraction)
    return k_values, mean_rho


def smooth_dense_k(
    k_coarse: np.ndarray,
    y_coarse: np.ndarray,
    n_dense: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Cubic spline resample for smooth panel curves."""
    if k_coarse.size < 4:
        k_dense = np.linspace(float(k_coarse[0]), float(k_coarse[-1]), n_dense)
        return k_dense, np.interp(k_dense, k_coarse, y_coarse)
    spl = CubicSpline(k_coarse, y_coarse, bc_type="natural")
    k_dense = np.linspace(float(k_coarse[0]), float(k_coarse[-1]), n_dense)
    return k_dense, spl(k_dense)


def mean_amplitude_trajectory(amps: np.ndarray) -> np.ndarray:
    r"""$\bar\rho(t)=N^{-1}\sum_j\rho_j(t)$."""
    return np.mean(amps, axis=0)


def make_figure(
    k_stable: np.ndarray,
    mean_rho_smooth: np.ndarray,
    k_theory: np.ndarray,
    m_eff_sq: np.ndarray,
    m_sq: float,
    t0: np.ndarray,
    rho_bar_lam0: np.ndarray,
    t1: np.ndarray,
    rho_bar_lam1: np.ndarray,
    *,
    prd_width: PrdWidth,
    burn_fraction: float,
    t_end: float,
    K_strong: float,
    lam_reg: float,
) -> plt.Figure:
    fig, (ax0, ax1, ax2) = plt.subplots(
        3,
        1,
        figsize=prd_figsize_three_panel(prd_width),
        layout="constrained",
    )

    ax0.plot(k_stable, mean_rho_smooth, color="C0")
    ax0.set_xlabel(r"Coupling $K$")
    ax0.set_ylabel(r"$\langle\rho\rangle$")
    ax0.set_xlim(0.0, float(k_stable[-1]))
    ax0.set_title(
        rf"(a) Mean amplitude, stable window $K<m^2$ (tail mean, burn ${100 * burn_fraction:.0f}\%$, "
        rf"$t_{{\mathrm{{f}}}}={t_end:g}$)"
    )
    ax0.grid(True, linestyle="--", alpha=0.55)

    ax1.plot(k_theory, m_eff_sq, color="C1")
    ax1.axhline(0.0, color="0.35", linestyle="--", linewidth=1.0, alpha=0.9)
    ax1.axvline(m_sq, color="0.35", linestyle=":", linewidth=1.2, alpha=0.9)
    ax1.axvspan(m_sq, float(k_theory[-1]), color="0.75", alpha=0.22, linewidth=0)
    ax1.set_xlabel(r"Coupling $K$")
    ax1.set_ylabel(r"$m_{\mathrm{eff}}^2$")
    ax1.set_title(
        r"(b) Linear theory $m_{\mathrm{eff}}^2=m^2-K$ ($\dot{\bar{\theta}}=0$); "
        r"$K>m^2$ shaded (tachyonic)"
    )
    ax1.grid(True, linestyle="--", alpha=0.55)
    ax1.set_xlim(float(k_theory[0]), float(k_theory[-1]))

    ax2.plot(t0, rho_bar_lam0, color="C3", linestyle="--", label=rf"$\lambda=0$ ($K={K_strong:g}$)")
    ax2.plot(t1, rho_bar_lam1, color="C0", linestyle="-", label=rf"$\lambda={lam_reg:g}$ ($K={K_strong:g}$)")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$\bar{\rho}(t)$")
    ax2.set_title(r"(c) Regularized vs.\ unregularized mean amplitude at strong coupling")
    ax2.legend(loc="best", framealpha=0.92)
    ax2.grid(True, linestyle="--", alpha=0.55)

    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_6",
        help="Output path without extension (default: code/figures/figure_6).",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument(
        "--prd-width",
        choices=("column", "full"),
        default="full",
        help="PRD canvas width: 'column' (~3.5 in) or 'full' (~7.0 in).",
    )
    p.add_argument("--t-end", type=float, default=45.0, help="Integration horizon for panel (a).")
    p.add_argument("--n-points", type=int, default=1200, help="Uniform time samples.")
    p.add_argument(
        "--burn-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction from the start before tail averaging (panel (a)).",
    )
    p.add_argument("--seed", type=int, default=42, help="RNG seed (match Figs.~1--3).")
    p.add_argument("--N", type=int, default=5, dest="N_fields", help="Number of oscillators.")
    p.add_argument("--m", type=float, default=1.0, help="Identical bare mass for each field.")
    p.add_argument("--k-stable-max", type=float, default=0.9, help="Upper K for panel (a).")
    p.add_argument("--n-k-stable", type=int, default=28, help="Coarse K samples in [0, k-stable-max].")
    p.add_argument("--n-k-smooth", type=int, default=220, help="Dense K samples for spline display.")
    p.add_argument("--k-theory-max", type=float, default=1.75, help="Upper K for panel (b).")
    p.add_argument("--n-k-theory", type=int, default=400, help="K samples for panel (b).")
    p.add_argument("--K-strong", type=float, default=8.0, help="Coupling for panel (c).")
    p.add_argument("--lam-reg", type=float, default=0.1, help="Quartic strength for regularized curve.")
    p.add_argument("--t-end-c", type=float, default=18.0, help="Time horizon for panel (c).")
    p.add_argument("--show", action="store_true", help="Show interactively.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    warnings.filterwarnings(
        "ignore",
        message="overflow encountered in square",
        category=RuntimeWarning,
    )
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    if not (0.0 <= args.burn_fraction < 1.0):
        raise SystemExit("error: --burn-fraction must satisfy 0 <= f < 1")
    if args.k_stable_max <= 0.0:
        raise SystemExit("error: --k-stable-max must be positive")
    if args.n_k_stable < 4:
        raise SystemExit("error: need at least 4 coarse K points for spline")

    N = args.N_fields
    masses = np.full(N, float(args.m), dtype=np.float64)
    m_sq = float(args.m) ** 2

    k_coarse = np.linspace(0.0, args.k_stable_max, args.n_k_stable)
    configure_matplotlib(prd_width)
    _k_c, mean_rho_c = scan_mean_rho_stable(
        k_coarse,
        N=N,
        masses=masses,
        t_end=args.t_end,
        n_points=args.n_points,
        seed=args.seed,
        burn_fraction=args.burn_fraction,
    )
    k_stable, mean_rho_smooth = smooth_dense_k(k_coarse, mean_rho_c, args.n_k_smooth)

    k_theory = np.linspace(0.0, args.k_theory_max, args.n_k_theory)
    m_eff_sq = m_sq - k_theory

    t0, amps0, ok0 = simulate(
        N,
        float(args.K_strong),
        masses,
        t_end=args.t_end_c,
        n_points=args.n_points,
        seed=args.seed,
        lam=0.0,
    )
    t1, amps1, ok1 = simulate(
        N,
        float(args.K_strong),
        masses,
        t_end=args.t_end_c,
        n_points=args.n_points,
        seed=args.seed,
        lam=float(args.lam_reg),
    )
    if not ok1:
        raise RuntimeError("panel (c): lambda>0 integration failed unexpectedly")

    rho0 = mean_amplitude_trajectory(amps0)
    rho1 = mean_amplitude_trajectory(amps1)
    if not ok0:
        warnings.warn("panel (c): lambda=0 run may be truncated near blow-up")

    fig = make_figure(
        k_stable,
        mean_rho_smooth,
        k_theory,
        m_eff_sq,
        m_sq,
        t0,
        rho0,
        t1,
        rho1,
        prd_width=prd_width,
        burn_fraction=args.burn_fraction,
        t_end=args.t_end,
        K_strong=float(args.K_strong),
        lam_reg=float(args.lam_reg),
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
