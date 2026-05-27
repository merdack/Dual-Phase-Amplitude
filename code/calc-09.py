r"""calc-09.py

Three-dimensional **finite-difference** toy integration of the coupled amplitude--phase
equations (14) and (16) from the manuscript, in the **Minkowski / local
patch** limit (\(a=1\), \(H=0\)) on a **periodic cubic** grid.

The code advances \(N\) scalar fields \(\{\rho_i,\theta_i\}_{i=1}^N\) with

* Eq.~(14) (real part): \(\ddot\rho_i - \nabla^2\rho_i - \rho_i(\dot\theta_i^2-|\nabla\theta_i|^2)
  + m_i^2\rho_i = (K/N)\sum_k \rho_k\cos(\theta_k-\theta_i)\).
* Eq.~(16) (imaginary part), solved for \(\ddot\theta_i\):
  \(\rho_i\ddot\theta_i = (K/N)\sum_k\rho_k\sin(\theta_k-\theta_i)
  -2(\dot\rho_i\dot\theta_i-\nabla\rho_i\!\cdot\!\nabla\theta_i)+\rho_i\nabla^2\theta_i\).

**Optional gravity (phenomenological, not derived from \(T_{\mu\nu}\) here):**

* ``external``: fixed soft Plummer \(\Phi(\mathbf{x})\); adds \(-\lambda_{\mathrm{grav}}\,\rho_i\,\Phi\)
  to the **amplitude** acceleration (toy ``mass in potential'' channel).
* ``self``: periodic Poisson \(\nabla^2\Phi = 4\pi G(\rho_{\mathrm{tot}}-\bar\rho_{\mathrm{tot}})\)
  via FFT each RHS evaluation, same coupling to \(\ddot\rho_i\).

This is a **research-grade toy** for qualitative comparison (e.g.\ \(K=0\) vs \(K>0\)),
not a cosmological N-body or GR solver.

Run::

    python calc-09.py [--gravity none|external|self] [--prd-width column|full] ...

Figure index matches script number: default ``figures/figure_9``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

PrdWidth = Literal["column", "full"]


def prd_figsize(width: PrdWidth) -> tuple[float, float]:
    if width == "column":
        return (7.2, 8.6)
    return (13.5, 9.0)


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


def laplacian_3d(f: np.ndarray, dx: float) -> np.ndarray:
    """7-point Laplacian on last three axes (periodic via ``np.roll``)."""
    inv_dx2 = 1.0 / (dx * dx)
    out = np.zeros_like(f, dtype=np.float64)
    for ax in (-3, -2, -1):
        out += (np.roll(f, 1, axis=ax) + np.roll(f, -1, axis=ax) - 2.0 * f) * inv_dx2
    return out


def grad_component(f: np.ndarray, axis: int, dx: float) -> np.ndarray:
    """Central gradient on ``axis`` (periodic)."""
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2.0 * dx)


def grad_squared_dot(
    gx_a: np.ndarray,
    gy_a: np.ndarray,
    gz_a: np.ndarray,
    gx_b: np.ndarray,
    gy_b: np.ndarray,
    gz_b: np.ndarray,
) -> np.ndarray:
    return gx_a * gx_b + gy_a * gy_b + gz_a * gz_b


def kuramoto_cos(rho: np.ndarray, theta: np.ndarray, K: float) -> np.ndarray:
    r"""Shape ``(N, nx, ny, nz)``: \((K/N)\sum_k \rho_k\cos(\theta_k-\theta_i)\)."""
    n = rho.shape[0]
    out = np.zeros_like(rho, dtype=np.float64)
    for i in range(n):
        for k in range(n):
            out[i] += rho[k] * np.cos(theta[k] - theta[i])
    return (K / float(n)) * out


def kuramoto_sin(rho: np.ndarray, theta: np.ndarray, K: float) -> np.ndarray:
    r"""\((K/N)\sum_k \rho_k\sin(\theta_k-\theta_i)\)."""
    n = rho.shape[0]
    out = np.zeros_like(rho, dtype=np.float64)
    for i in range(n):
        for k in range(n):
            out[i] += rho[k] * np.sin(theta[k] - theta[i])
    return (K / float(n)) * out


def total_density_modulus_squared(rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
    r"""\(|\sum_i \rho_i e^{i\theta_i}|^2\) on the grid."""
    z = np.zeros(rho.shape[1:], dtype=np.complex128)
    for i in range(rho.shape[0]):
        z += rho[i] * np.exp(1j * theta[i])
    return np.abs(z) ** 2


def external_plummer_potential(shape: tuple[int, int, int], dx: float, G: float, M: float, eps: float) -> np.ndarray:
    nx, ny, nz = shape
    ix = np.arange(nx, dtype=np.float64) - 0.5 * (nx - 1)
    iy = np.arange(ny, dtype=np.float64) - 0.5 * (ny - 1)
    iz = np.arange(nz, dtype=np.float64) - 0.5 * (nz - 1)
    X, Y, Z = np.meshgrid(ix * dx, iy * dx, iz * dx, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z + eps * eps)
    return (-G * M / r).astype(np.float64)


def poisson_periodic_fft(
    source: np.ndarray,
    dx: float,
    G_newton: float,
) -> np.ndarray:
    r"""Solve \(\nabla^2\Phi = 4\pi G\,S\) with periodic BC (FFT), \(\hat\Phi(\mathbf k=0)=0\)."""
    nx, ny, nz = source.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dx)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    k2 = KX * KX + KY * KY + KZ * KZ
    src_hat = np.fft.fftn(source)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi_hat = np.where(k2 > 0.0, 4.0 * np.pi * G_newton * src_hat / k2, 0.0)
    phi_hat[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(phi_hat)).astype(np.float64)


def time_derivatives(
    rho: np.ndarray,
    pi_rho: np.ndarray,
    theta: np.ndarray,
    pi_theta: np.ndarray,
    *,
    K: float,
    m: np.ndarray,
    dx: float,
    Phi: np.ndarray | None,
    coupl_grav: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Return \((\dot\rho,\ddot\rho,\dot\theta,\ddot\theta)\) arrays."""
    lap_rho = laplacian_3d(rho, dx)
    lap_theta = laplacian_3d(theta, dx)

    gx_r = grad_component(rho, -3, dx)
    gy_r = grad_component(rho, -2, dx)
    gz_r = grad_component(rho, -1, dx)
    gx_t = grad_component(theta, -3, dx)
    gy_t = grad_component(theta, -2, dx)
    gz_t = grad_component(theta, -1, dx)

    grad_sq = grad_squared_dot(gx_t, gy_t, gz_t, gx_t, gy_t, gz_t)
    cross = grad_squared_dot(gx_r, gy_r, gz_r, gx_t, gy_t, gz_t)

    kcos = kuramoto_cos(rho, theta, K)
    ksin = kuramoto_sin(rho, theta, K)

    ddot_rho = (
        lap_rho
        - rho * (pi_theta * pi_theta - grad_sq)
        + (m[:, None, None, None] ** 2) * rho
        + kcos
    )
    if Phi is not None and coupl_grav != 0.0:
        ddot_rho -= coupl_grav * rho * Phi

    rho_safe = np.maximum(rho, 1e-8)
    A = pi_rho * pi_theta - cross
    ddot_theta = (ksin - 2.0 * A + rho * lap_theta) / rho_safe

    return pi_rho, ddot_rho, pi_theta, ddot_theta


def pack_state(rho: np.ndarray, pi_rho: np.ndarray, theta: np.ndarray, pi_theta: np.ndarray) -> np.ndarray:
    return np.concatenate((rho.ravel(), pi_rho.ravel(), theta.ravel(), pi_theta.ravel()))


def unpack_state(
    y: np.ndarray,
    *,
    n: int,
    nx: int,
    ny: int,
    nz: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    size = n * nx * ny * nz
    rho = y[0:size].reshape(n, nx, ny, nz)
    pi_rho = y[size : 2 * size].reshape(n, nx, ny, nz)
    theta = y[2 * size : 3 * size].reshape(n, nx, ny, nz)
    pi_theta = y[3 * size :].reshape(n, nx, ny, nz)
    return rho, pi_rho, theta, pi_theta


def rk4_step_fixed(
    y: np.ndarray,
    dt: float,
    *,
    K: float,
    m: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
    dx: float,
    Phi_ext: np.ndarray | None,
    G_newton: float,
    coupl_grav: float,
    gravity: str,
) -> np.ndarray:
    n = m.size

    def rhs(state: np.ndarray) -> np.ndarray:
        rho, pr, th, pt = unpack_state(state, n=n, nx=nx, ny=ny, nz=nz)
        phi: np.ndarray | None = Phi_ext
        if gravity == "self":
            rho_tot = np.sum(rho, axis=0)
            src = rho_tot - float(np.mean(rho_tot))
            phi = poisson_periodic_fft(src, dx, G_newton)
        dr, ddr, dt1, ddt = time_derivatives(rho, pr, th, pt, K=K, m=m, dx=dx, Phi=phi, coupl_grav=coupl_grav)
        return pack_state(dr, ddr, dt1, ddt)

    k1 = rhs(y)
    k2 = rhs(y + 0.5 * dt * k1)
    k3 = rhs(y + 0.5 * dt * k2)
    k4 = rhs(y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def initial_data(
    *,
    n: int,
    nx: int,
    ny: int,
    nz: int,
    dx: float,
    rng: np.random.Generator,
    rho_bg: float,
    bump: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rho = np.full((n, nx, ny, nz), rho_bg, dtype=np.float64)
    theta = rng.uniform(-0.15, 0.15, size=(n, nx, ny, nz)).astype(np.float64)
    ix = np.arange(nx, dtype=np.float64) - 0.5 * (nx - 1)
    iy = np.arange(ny, dtype=np.float64) - 0.5 * (ny - 1)
    iz = np.arange(nz, dtype=np.float64) - 0.5 * (nz - 1)
    X, Y, Z = np.meshgrid(ix * dx, iy * dx, iz * dx, indexing="ij")
    for i in range(n):
        phase_shift = 2 * np.pi * i / max(n, 1)
        r2 = (X - 0.8 * dx * np.cos(phase_shift)) ** 2 + (Y - 0.6 * dx * np.sin(phase_shift)) ** 2 + Z * Z
        rho[i] += bump * np.exp(-r2 / (2.5 * dx) ** 2)
    pi_rho = rng.normal(0.0, 0.02 * bump, size=(n, nx, ny, nz)).astype(np.float64)
    pi_theta = rng.normal(0.0, 0.02, size=(n, nx, ny, nz)).astype(np.float64)
    return rho, pi_rho, theta, pi_theta


def run_simulation(
    *,
    K: float,
    n: int,
    nx: int,
    ny: int,
    nz: int,
    dx: float,
    dt: float,
    n_steps: int,
    m: np.ndarray,
    rng: np.random.Generator,
    gravity: str,
    G_newton: float,
    M_ext: float,
    eps_plummer: float,
    coupl_grav: float,
    rho_bg: float,
    bump: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(rho_mod2_history, max_mod2)`` shapes ``(n_steps+1, nx, ny, nz)`` and ``(n_steps+1,)``."""
    rho, pi_rho, theta, pi_theta = initial_data(
        n=n, nx=nx, ny=ny, nz=nz, dx=dx, rng=rng, rho_bg=rho_bg, bump=bump
    )
    Phi_ext: np.ndarray | None = None
    if gravity == "external":
        Phi_ext = external_plummer_potential((nx, ny, nz), dx, G_newton, M_ext, eps_plummer)
    elif gravity == "none":
        Phi_ext = None
    elif gravity == "self":
        Phi_ext = None
    else:
        raise ValueError("gravity must be none|external|self")

    y = pack_state(rho, pi_rho, theta, pi_theta)
    hist = np.empty((n_steps + 1, nx, ny, nz), dtype=np.float64)
    mx = np.empty(n_steps + 1, dtype=np.float64)
    rho0, pr0, th0, pt0 = unpack_state(y, n=n, nx=nx, ny=ny, nz=nz)
    hist[0] = total_density_modulus_squared(rho0, th0)
    mx[0] = float(np.max(hist[0]))
    for s in range(n_steps):
        y = rk4_step_fixed(
            y,
            dt,
            K=K,
            m=m,
            nx=nx,
            ny=ny,
            nz=nz,
            dx=dx,
            Phi_ext=Phi_ext,
            G_newton=G_newton,
            coupl_grav=coupl_grav,
            gravity=gravity,
        )
        rho, pr, th, pt = unpack_state(y, n=n, nx=nx, ny=ny, nz=nz)
        hist[s + 1] = total_density_modulus_squared(rho, th)
        mx[s + 1] = float(np.max(hist[s + 1]))
    return hist, mx


def make_figure(
    hist0: np.ndarray,
    histk: np.ndarray,
    mx0: np.ndarray,
    mxk: np.ndarray,
    t: np.ndarray,
    z_index: int,
    *,
    K0: float,
    Kk: float,
    gravity: str,
    prd_width: PrdWidth,
) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=prd_figsize(prd_width))
    sl0 = hist0[-1, :, :, z_index]
    slk = histk[-1, :, :, z_index]
    vmax = float(np.percentile(np.concatenate((sl0.ravel(), slk.ravel())), 99.5))
    vmax = max(vmax, 1e-18)

    im0 = axes[0, 0].imshow(
        sl0.T,
        origin="lower",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        aspect="auto",
    )
    axes[0, 0].set_title(rf"(A) $|\Phi|^2$ slice ($K={K0:g}$, gravity={gravity})")
    axes[0, 0].set_xlabel(r"$x$ index")
    axes[0, 0].set_ylabel(r"$y$ index")
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)

    im1 = axes[0, 1].imshow(
        slk.T,
        origin="lower",
        cmap="magma",
        vmin=0.0,
        vmax=vmax,
        aspect="auto",
    )
    axes[0, 1].set_title(rf"(B) $|\Phi|^2$ slice ($K={Kk:g}$, gravity={gravity})")
    axes[0, 1].set_xlabel(r"$x$ index")
    axes[0, 1].set_ylabel(r"$y$ index")
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04)

    axes[1, 0].plot(t, mx0, label=rf"$K={K0:g}$")
    axes[1, 0].plot(t, mxk, label=rf"$K={Kk:g}$")
    axes[1, 0].set_xlabel(r"time $t$")
    axes[1, 0].set_ylabel(r"$\max |\Phi|^2$")
    axes[1, 0].set_title(r"(C) Peak interference density vs time")
    axes[1, 0].grid(True, linestyle="--", alpha=0.55)
    axes[1, 0].legend(loc="best")

    v0 = float(np.var(mx0[10:])) if mx0.size > 15 else float(np.var(mx0))
    vk = float(np.var(mxk[10:])) if mxk.size > 15 else float(np.var(mxk))
    axes[1, 1].axis("off")
    axes[1, 1].text(
        0.02,
        0.75,
        rf"Temporal Var($\max|\Phi|^2$) after step 10:" "\n"
        rf"$K={K0:g}$: {v0:.4g}" "\n"
        rf"$K={Kk:g}$: {vk:.4g}" "\n\n"
        r"Discretized Eqs.~(14)--(16), $a=1$, $H=0$, periodic cube." "\n"
        r"Optional gravity couples as $-\lambda \rho_i \Phi$ in $\ddot\rho_i$ (see docstring).",
        fontsize=10,
        va="top",
        family="monospace",
    )

    fig.subplots_adjust(hspace=0.35, wspace=0.35, top=0.93, bottom=0.07)
    return fig


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="3D finite-difference demo of coupled Eqs. (14)+(16) with optional toy gravity.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures" / "figure_9",
        help="Output basename without extension.",
    )
    p.add_argument("--format", default="pdf,png", help="Comma-separated formats.")
    p.add_argument("--prd-width", choices=("column", "full"), default="full")
    p.add_argument("--N", type=int, default=4, help="Number of coupled scalar fields.")
    p.add_argument("--nx", type=int, default=20, help="Grid size (cubic: nx=ny=nz).")
    p.add_argument("--dt", type=float, default=0.012, help="RK4 time step.")
    p.add_argument("--n-steps", type=int, default=45, help="Number of RK4 steps.")
    p.add_argument(
        "--dx",
        type=float,
        default=0.35,
        help="Cell spacing; periodic cube linear extent is approximately nx*dx.",
    )
    p.add_argument("--K0", type=float, default=0.0, help="Coupling for competitor run.")
    p.add_argument("--K-kkg", type=float, default=3.5, dest="K_kkg", help="Coupling for KKG run.")
    p.add_argument("--m", type=float, default=1.0, help="Mass parameter (same for each field unless --m-spread).")
    p.add_argument(
        "--m-spread",
        type=float,
        default=0.04,
        help="Relative spread of masses across fields: m_i = m*(1+spread*(i/(N-1)-0.5)).",
    )
    p.add_argument("--seed", type=int, default=4242, help="RNG seed.")
    p.add_argument("--rho-bg", type=float, default=0.35, dest="rho_bg")
    p.add_argument("--bump", type=float, default=0.25, help="Gaussian bump amplitude on rho_i.")
    p.add_argument(
        "--gravity",
        choices=("none", "external", "self"),
        default="none",
        help="none | external Plummer Phi | self periodic Poisson Phi.",
    )
    p.add_argument("--G-newton", type=float, default=0.08, dest="G_newton", help="Newton constant for Poisson / Plummer.")
    p.add_argument("--M-ext", type=float, default=2.5, dest="M_ext", help="External Plummer mass parameter.")
    p.add_argument("--eps-plummer", type=float, default=1.2, dest="eps_plummer", help="Softening length (grid units * dx).")
    p.add_argument(
        "--coupl-grav",
        type=float,
        default=0.35,
        dest="coupl_grav",
        help="Strength of -lambda rho_i Phi term in ddot(rho).",
    )
    p.add_argument("--show", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prd_width: PrdWidth = "column" if args.prd_width == "column" else "full"

    n = int(args.N)
    nx = ny = nz = int(args.nx)
    if n < 2 or nx < 8:
        raise SystemExit("error: require N>=2 and nx>=8")

    rng = np.random.default_rng(args.seed)
    if n > 1:
        m = args.m * (1.0 + args.m_spread * (np.arange(n, dtype=np.float64) / (n - 1) - 0.5))
    else:
        m = np.array([args.m], dtype=np.float64)

    hist0, mx0 = run_simulation(
        K=args.K0,
        n=n,
        nx=nx,
        ny=ny,
        nz=nz,
        dx=args.dx,
        dt=args.dt,
        n_steps=args.n_steps,
        m=m,
        rng=rng,
        gravity=args.gravity,
        G_newton=args.G_newton,
        M_ext=args.M_ext,
        eps_plummer=args.eps_plummer,
        coupl_grav=args.coupl_grav,
        rho_bg=args.rho_bg,
        bump=args.bump,
    )

    rng2 = np.random.default_rng(args.seed)
    histk, mxk = run_simulation(
        K=args.K_kkg,
        n=n,
        nx=nx,
        ny=ny,
        nz=nz,
        dx=args.dx,
        dt=args.dt,
        n_steps=args.n_steps,
        m=m,
        rng=rng2,
        gravity=args.gravity,
        G_newton=args.G_newton,
        M_ext=args.M_ext,
        eps_plummer=args.eps_plummer,
        coupl_grav=args.coupl_grav,
        rho_bg=args.rho_bg,
        bump=args.bump,
    )

    t = np.linspace(0.0, args.dt * args.n_steps, args.n_steps + 1)
    z_index = nz // 2

    configure_matplotlib(prd_width)
    fig = make_figure(
        hist0,
        histk,
        mx0,
        mxk,
        t,
        z_index,
        K0=args.K0,
        Kk=args.K_kkg,
        gravity=args.gravity,
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
