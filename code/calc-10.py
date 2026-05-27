r"""calc-10.py

Toy Tolman--Oppenheimer--Volkoff (TOV) integration in **geometric units** ($G=c=1$) to
compare a **baseline barotropic** equation of state (stand-in for a self-interacting /
polytropic scalar condensate) against a **KKG-motivated stiffening** term
$P_{\mathrm{sync}}\propto K\rho^2$ added to the same baseline.

Parameters use **geometric units** ($G=c=1$); the baseline polytrope coefficients fix a single mass scale in this schematic (analogous to setting the boson mass $m$ to unity in a dimensionless closure).
This is **not** a microphysical derivation from the covariant KKG stress tensor; it is a
schematic barotropic closure $P(\rho;K)$ designed so that increasing Kuramoto-like
coupling $K$ raises the effective stiffness and shifts the mass--radius envelope toward
higher maximum masses in the same TOV channel.

The TOV system (Schwarzschild metric, static spherical symmetry) is
$\mathrm{d}M/\mathrm{d}r = 4\pi r^2\rho$ and
$\mathrm{d}P/\mathrm{d}r = -(\rho+P)(M+4\pi r^3 P)/(r(r-2M))$.

Run::

    python calc-10.py [--prd-width column|full] [--out PATH] [--format pdf,png] [--show]

Figure index matches script number: default ``figures/figure_10``.
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


def prd_figsize(width: PrdWidth) -> tuple[float, float]:
    if width == "column":
        return (3.5, 6.2)
    return (7.2, 6.2)


def configure_matplotlib(prd_width: PrdWidth) -> None:
    label = 11 if prd_width == "column" else 14
    tick = 10 if prd_width == "column" else 12
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
        }
    )


def pressure_eos(
    rho: float | np.ndarray,
    *,
    K: float,
    kappa: float,
    beta: float,
    alpha_sync: float,
) -> float | np.ndarray:
    r"""Barotropic pressure $P = \kappa \rho^{\beta} + \alpha_{\mathrm{sync}} K \rho^2$."""
    return kappa * np.power(rho, beta) + alpha_sync * K * np.power(rho, 2.0)


def rho_from_pressure(
    P_target: float,
    *,
    K: float,
    kappa: float,
    beta: float,
    alpha_sync: float,
    rho_hint: float | None = None,
) -> float:
    """Invert monotone $P(\rho)$ by Newton iterations (fast inside RK steps)."""
    if P_target <= 0.0:
        return 0.0

    a = alpha_sync * K
    rho = rho_hint if rho_hint is not None else (P_target / kappa) ** (1.0 / beta)
    rho = float(np.clip(rho, 1e-16, 1e8))

    for _ in range(48):
        Pv = float(kappa * rho**beta + a * rho**2)
        dP = float(kappa * beta * rho ** (beta - 1.0) + 2.0 * a * rho)
        if dP <= 0.0 or not np.isfinite(dP):
            break
        step = (Pv - P_target) / dP
        rho_new = rho - step
        if not np.isfinite(rho_new):
            break
        if abs(step) < 1e-14 * max(rho, 1.0):
            return float(max(rho_new, 1e-16))
        rho = float(np.clip(rho_new, 1e-16, 1e8))

    return float(max(rho, 1e-16))


def tov_rhs(
    r: float,
    y: np.ndarray,
    *,
    K: float,
    kappa: float,
    beta: float,
    alpha_sync: float,
) -> np.ndarray:
    """Derivatives for state ``y = [M, P]``."""
    M, P = float(y[0]), float(y[1])
    rho = rho_from_pressure(P, K=K, kappa=kappa, beta=beta, alpha_sync=alpha_sync)
    num = (rho + P) * (M + 4.0 * np.pi * r**3 * P)
    den = r * (r - 2.0 * M)
    if den <= 1e-14 * r**2 or r <= 0.0:
        return np.array([0.0, 0.0], dtype=np.float64)
    dM_dr = 4.0 * np.pi * r**2 * rho
    dP_dr = -num / den
    return np.array([dM_dr, dP_dr], dtype=np.float64)


def integrate_star(
    rho_c: float,
    *,
    K: float,
    kappa: float,
    beta: float,
    alpha_sync: float,
    r_min: float,
    r_max: float,
    P_surface: float,
) -> tuple[float, float, bool]:
    """Return ``(M_km, R_km, ok)`` in geometric length units (same as $M$ here)."""

    P_c = float(pressure_eos(rho_c, K=K, kappa=kappa, beta=beta, alpha_sync=alpha_sync))
    M0 = (4.0 * np.pi / 3.0) * r_min**3 * rho_c
    y0 = np.array([M0, P_c], dtype=np.float64)

    def rhs(r: float, y: np.ndarray) -> np.ndarray:
        return tov_rhs(r, y, K=K, kappa=kappa, beta=beta, alpha_sync=alpha_sync)

    def reached_surface(r: float, y: np.ndarray) -> float:
        return float(y[1]) - P_surface

    reached_surface.terminal = True
    reached_surface.direction = -1.0

    def near_horizon(r: float, y: np.ndarray) -> float:
        M, _P = float(y[0]), float(y[1])
        return 0.99 * r - 2.0 * M

    near_horizon.terminal = True
    near_horizon.direction = -1.0

    sol = solve_ivp(
        rhs,
        (r_min, r_max),
        y0,
        method="RK45",
        rtol=1e-7,
        atol=1e-9,
        events=[reached_surface, near_horizon],
        dense_output=False,
        max_step=r_max * 0.02,
    )

    if not sol.success or sol.t.size < 2:
        return 0.0, 0.0, False

    # Surface: last event where P hit P_surface, else last point with P>P_surface
    R = float(sol.t[-1])
    M = float(sol.y[0, -1])
    if sol.t_events[0] is not None and sol.t_events[0].size > 0:
        R = float(sol.t_events[0][-1])
        M = float(sol.y_events[0][-1, 0])
    if sol.t_events[1] is not None and sol.t_events[1].size > 0:
        return M, R, False
    if M <= 0.0 or R <= r_min:
        return M, R, False
    return M, R, True


def mass_radius_sequence(
    rho_cs: np.ndarray,
    *,
    K: float,
    kappa: float,
    beta: float,
    alpha_sync: float,
    r_min: float,
    r_max: float,
    P_surface: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(R, M, ok)`` arrays for each central density."""
    R_list: list[float] = []
    M_list: list[float] = []
    ok_list: list[bool] = []
    for rho_c in rho_cs:
        M, R, ok = integrate_star(
            float(rho_c),
            K=K,
            kappa=kappa,
            beta=beta,
            alpha_sync=alpha_sync,
            r_min=r_min,
            r_max=r_max,
            P_surface=P_surface,
        )
        R_list.append(R)
        M_list.append(M)
        ok_list.append(ok)
    return (
        np.asarray(R_list, dtype=np.float64),
        np.asarray(M_list, dtype=np.float64),
        np.asarray(ok_list, dtype=bool),
    )


def stable_branch_mask(M: np.ndarray, ok: np.ndarray) -> np.ndarray:
    """Keep successful integrations up to the first local maximum of ``M`` versus central density."""
    out = np.zeros_like(ok, dtype=bool)
    valid = np.where(ok & np.isfinite(M) & (M > 0.0))[0]
    if valid.size == 0:
        return out
    Mv = M[valid]
    peak_local = Mv.size - 1
    for t in range(1, Mv.size):
        if Mv[t] < Mv[t - 1] - 1e-10:
            peak_local = t - 1
            break
    kept = valid[: peak_local + 1]
    out[kept] = True
    return out


def m_max_radius_at_max(
    R: np.ndarray,
    M: np.ndarray,
    mask: np.ndarray,
) -> tuple[float, float]:
    if not np.any(mask):
        return 0.0, 0.0
    i = int(np.argmax(M[mask]))
    w = np.where(mask)[0][i]
    return float(M[w]), float(R[w])


def make_figure(
    curves: list[tuple[float, np.ndarray, np.ndarray, np.ndarray]],
    mmax_K: list[tuple[float, float]],
    *,
    prd_width: PrdWidth,
) -> plt.Figure:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=prd_figsize(prd_width))

    colors = plt.cm.viridis(np.linspace(0.15, 0.95, max(len(curves), 2)))
    for c, (K, R, M, mask) in enumerate(curves):
        ax0.plot(R[mask], M[mask], color=colors[c], lw=2.0, label=rf"$K={K:g}$")
        mm, rm = m_max_radius_at_max(R, M, mask)
        if mm > 0.0:
            ax0.scatter([rm], [mm], color=colors[c], s=38, zorder=5, edgecolors="k", linewidths=0.4)

    ax0.set_xlabel(r"radius $R$ ($G=c=1$)")
    ax0.set_ylabel(r"mass $M$ ($G=c=1$)")
    ax0.set_title(r"(A) Mass--radius sequences (toy EoS)")
    ax0.grid(True, linestyle="--", alpha=0.45)
    ax0.legend(loc="best", fontsize=8)

    Ks = np.array([k for k, _ in mmax_K], dtype=np.float64)
    Mms = np.array([mm for _, mm in mmax_K], dtype=np.float64)
    ax1.plot(Ks, Mms, "o-", color="darkred", lw=2.0, markersize=6)
    ax1.set_xlabel(r"Kuramoto-like coupling $K$")
    ax1.set_ylabel(r"$M_{\max}$ ($G=c=1$)")
    ax1.set_title(r"(B) Peak TOV mass vs $K$")
    ax1.grid(True, linestyle="--", alpha=0.45)

    fig.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.12, wspace=0.32)
    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Toy TOV mass--radius curves: baseline EoS vs K-stiffened EoS.")
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_10",
        help="Output basename without extension.",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument("--prd-width", choices=("column", "full"), default="full")
    p.add_argument("--seed", type=int, default=7, help="Reserved for reproducibility hooks.")
    p.add_argument(
        "--K-list",
        type=str,
        default="0,2,4,6,8,10",
        help="Comma-separated K values for M-R curves.",
    )
    p.add_argument("--rho-min", type=float, default=0.002, dest="rho_min", help="Min central density (code units).")
    p.add_argument("--rho-max", type=float, default=0.14, dest="rho_max", help="Max central density (code units).")
    p.add_argument("--n-rho", type=int, default=28, dest="n_rho", help="Number of central densities (log-spaced).")
    p.add_argument("--kappa", type=float, default=12.0, help=r"Polytrope prefactor $\kappa$ in $P=\kappa\rho^\beta+\cdots$.")
    p.add_argument("--beta", type=float, default=1.8, help=r"Polytrope exponent $\beta$ (must exceed unity).")
    p.add_argument(
        "--alpha-sync",
        type=float,
        default=2.2,
        dest="alpha_sync",
        help=r"Coefficient $\alpha_{\mathrm{sync}}$ in $P_{\mathrm{sync}}=\alpha_{\mathrm{sync}} K \rho^2$.",
    )
    p.add_argument("--r-min", type=float, default=5e-5, dest="r_min", help="Inner starting radius.")
    p.add_argument("--r-max", type=float, default=80.0, dest="r_max", help="Fallback outer radius cap.")
    p.add_argument("--P-surface", type=float, default=1e-10, dest="P_surface", help="Pressure cutoff at surface.")
    p.add_argument("--show", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"
    _ = np.random.default_rng(args.seed)

    K_list = [float(x.strip()) for x in args.K_list.split(",") if x.strip()]
    rho_cs = np.geomspace(args.rho_min, args.rho_max, args.n_rho)

    curves: list[tuple[float, np.ndarray, np.ndarray, np.ndarray]] = []
    mmax_K: list[tuple[float, float]] = []

    for K in K_list:
        R, M, ok = mass_radius_sequence(
            rho_cs,
            K=K,
            kappa=args.kappa,
            beta=args.beta,
            alpha_sync=args.alpha_sync,
            r_min=args.r_min,
            r_max=args.r_max,
            P_surface=args.P_surface,
        )
        mask = stable_branch_mask(M, ok)
        curves.append((K, R, M, mask))
        mm, _rm = m_max_radius_at_max(R, M, mask)
        mmax_K.append((K, mm))

    configure_matplotlib(prd_width)
    fig = make_figure(curves, mmax_K, prd_width=prd_width)

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
