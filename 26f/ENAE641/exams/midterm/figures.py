"""Figures for the ENAE 641 midterm formula sheet.

Writes the conceptual plots (SVG) and the lecture-slide crops (PNG) into ./assets:

    nix develop --command python figures.py
"""

import subprocess
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

matplotlib.use("svg")

ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "assets"
LECTURES = ROOT.parents[1] / "docs/lectures"

# categorical slots 1-3 (validated all-pairs), status colours (always paired with a text label), chart ink
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GOOD, SERIOUS, CRITICAL = "#0ca30c", "#ec835a", "#d03b3b"
INK, SECONDARY, MUTED = "#0b0b0b", "#52514e", "#898781"
AXIS = "#c3c2b7"
plt.rcParams.update(
    {
        "axes.edgecolor": AXIS,
        "axes.formatter.use_mathtext": True,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.5,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.titlesize": 6.5,
        "axes.titlepad": 2,
        "font.family": "serif",
        "font.serif": ["cmr10"],
        "font.size": 6.5,
        "lines.linewidth": 0.9,
        "lines.solid_capstyle": "round",
        "mathtext.fontset": "cm",
        "savefig.pad_inches": 0.01,
        "svg.fonttype": "path",
        "text.color": INK,
        "xtick.color": MUTED,
        "xtick.labelsize": 5.5,
        "xtick.major.size": 1.5,
        "xtick.major.width": 0.4,
        "ytick.color": MUTED,
        "ytick.labelsize": 5.5,
        "ytick.major.size": 1.5,
        "ytick.major.width": 0.4,
    }
)


def save(fig, name):
    fig.savefig(ASSETS / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def bare(ax):
    """Strip ticks and spines (schematic axes)."""
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def flow(A, x0, t):
    """x(t) = e^{At} x0 via the eigendecomposition (distinct eigenvalues)."""
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.asarray(x0, dtype=complex))
    return (V @ (c[:, None] * np.exp(np.outer(w, t)))).real


def arrow(ax, start, delta, color, lw=0.7):
    ax.annotate(
        "",
        xy=start + delta,
        xytext=start,
        arrowprops={"arrowstyle": "-|>", "color": color, "lw": lw, "mutation_scale": 5, "shrinkA": 0, "shrinkB": 0},
    )


def splane_modes():
    """Pole location in the s-plane vs. the shape of its (real) mode e^{σt} cos ωt."""
    fig, ax = plt.subplots(figsize=(2.5, 1.5))
    ax.set_xlim(-3.3, 2.15)
    ax.set_ylim(-0.62, 3.25)
    bare(ax)
    ax.axvspan(-3.3, 0, color=GOOD, alpha=0.08, lw=0)
    ax.axhline(0, color=AXIS, lw=0.5, zorder=0)
    ax.axvline(0, color=AXIS, lw=0.5, zorder=0)
    ax.text(2.13, -0.06, r"$\sigma=\mathrm{Re}\,p$", ha="right", va="top", color=SECONDARY, fontsize=5.5)
    ax.text(0.06, 3.22, r"$\omega=\mathrm{Im}\,p$", ha="left", va="top", color=SECONDARY, fontsize=5.5)

    # each inset shows the mode's shape only (normalized to its largest value on 0 <= t <= 6)
    t = np.linspace(0, 6, 400)
    for sigma, color in [(-2.4, GOOD), (-1.0, GOOD), (0.0, SERIOUS), (1.15, CRITICAL)]:
        for omega, y0 in [(0.0, 0.14), (1.8, 1.95)]:
            ax.plot(sigma, omega, marker="x", ms=3.2, mew=0.9, color=color, zorder=3)
            ins = ax.inset_axes([sigma - 0.47, y0, 0.94, 0.78], transform=ax.transData)
            mode = np.exp(0.35 * sigma * t) * np.cos(omega * t)
            ins.plot(t, mode / np.abs(mode).max(), color=color, lw=0.8)
            ins.axhline(0, color=AXIS, lw=0.4)
            ins.set_xlim(0, 6)
            ins.set_ylim(-1.15, 1.15)
            bare(ins)
            ins.patch.set_alpha(0)
            for spine in ins.spines.values():
                spine.set_visible(True)
                spine.set_color(AXIS)
                spine.set_linewidth(0.3)

    for x, label in [(-1.7, "stable (LHP)"), (0.0, "marginal"), (1.15, "unstable")]:
        ax.text(x, -0.47, label, ha="center", va="center", color=SECONDARY, fontsize=5.5)
    save(fig, "splane-modes")


def trdet():
    """Trace-determinant plane: phase-portrait type of a 2x2 LTI system."""
    fig, ax = plt.subplots(figsize=(2.5, 1.8))
    ax.set_xlim(-4.3, 4.3)
    ax.set_ylim(-2.75, 5.55)
    bare(ax)
    tau = np.linspace(-4.3, 4.3, 400)
    ax.fill_between(tau, 0, 5.55, where=tau <= 0, color=GOOD, alpha=0.08, lw=0)
    ax.fill_between(tau, 0, 5.55, where=tau >= 0, color=CRITICAL, alpha=0.06, lw=0)
    ax.fill_between(tau, -2.75, 0, color=CRITICAL, alpha=0.06, lw=0)
    ax.axhline(0, color=AXIS, lw=0.5, zorder=0)
    ax.axvline(0, color=AXIS, lw=0.5, zorder=0)
    ax.plot([0, 0], [0, 5.55], color=SERIOUS, lw=1.2, zorder=1)
    ax.plot(tau, tau**2 / 4, color=SECONDARY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.text(4.25, -0.12, r"$\tau$", ha="right", va="top", color=SECONDARY, fontsize=6)
    ax.text(0.12, 2.3, r"$\Delta$", ha="left", va="center", color=SECONDARY, fontsize=6)
    ax.text(-4.2, 3.95, r"$\Delta=\tau^2/4$", ha="left", va="center", color=SECONDARY, fontsize=5.5)

    g = np.linspace(-1, 1, 25)
    X, Y = np.meshgrid(g, g)
    cases = [
        ("stable node", [[-2.0, 0.6], [0.3, -1.0]], (-3.3, 1.47), "below"),
        ("stable spiral", [[-0.5, 2.0], [-2.0, -0.5]], (-1.9, 4.3), "below"),
        ("center", [[0.0, 2.0], [-1.2, 0.0]], (0.0, 4.3), "below"),
        ("unstable spiral", [[0.5, 2.0], [-2.0, 0.5]], (1.9, 4.3), "below"),
        ("unstable node", [[2.0, 0.6], [0.3, 1.0]], (3.3, 1.47), "below"),
        ("saddle", [[1.0, 0.5], [0.8, -1.0]], (0.0, -1.6), "right"),
    ]
    for name, A, (cx, cy), where in cases:
        A = np.array(A)
        ins = ax.inset_axes([cx - 0.72, cy - 1.02, 1.44, 2.04], transform=ax.transData)
        ins.streamplot(
            X, Y, A[0, 0] * X + A[0, 1] * Y, A[1, 0] * X + A[1, 1] * Y,
            density=0.42, color=INK, linewidth=0.35, arrowsize=0.35,
        )
        w, V = np.linalg.eig(A)
        if np.isreal(w).all():
            for v in V.T:
                v = v.real / np.abs(v.real).max()
                ins.plot([-v[0], v[0]], [-v[1], v[1]], color=BLUE, lw=0.6)
        ins.set_xlim(-1, 1)
        ins.set_ylim(-1, 1)
        bare(ins)
        ins.set_facecolor("white")
        for spine in ins.spines.values():
            spine.set_visible(True)
            spine.set_color(AXIS)
            spine.set_linewidth(0.3)
        if where == "below":
            ax.text(cx, cy - 1.1, name, ha="center", va="top", fontsize=5.5)
        else:
            ax.text(cx + 0.82, cy, name, ha="left", va="center", fontsize=5.5)
    save(fig, "trdet")


def pendulum_cases():
    """Pendulum (g/l = 5, b'/(ml^2) = 2) released 0.3 rad short of two equilibria."""
    fig, axs = plt.subplots(1, 2, figsize=(2.5, 1.0))
    t = np.linspace(0, 7.5, 600)
    y0 = -0.3

    # θ* = π/4: poles -1 ± j ω_d, free response from y(0) = y0, y'(0) = 0
    wd = np.sqrt(5 * np.cos(np.pi / 4) - 1)
    y = y0 * np.exp(-t) * (np.cos(wd * t) + np.sin(wd * t) / wd)
    env = abs(y0) * np.sqrt(1 + 1 / wd**2) * np.exp(-t)
    ax = axs[0]
    ax.axhline(45, color=MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.plot(t, 45 + np.degrees(env), color=SECONDARY, lw=0.5, ls=(0, (1, 1.5)))
    ax.plot(t, 45 - np.degrees(env), color=SECONDARY, lw=0.5, ls=(0, (1, 1.5)))
    ax.plot(t, 45 + np.degrees(y), color=BLUE)
    t1, t2 = np.pi / wd, 3 * np.pi / wd
    ax.annotate("", xy=(t2, 55), xytext=(t1, 55), arrowprops={"arrowstyle": "<->", "lw": 0.5, "color": SECONDARY, "mutation_scale": 4, "shrinkA": 0, "shrinkB": 0})
    ax.text((t1 + t2) / 2, 55.8, r"$2\pi/\omega_d$", ha="center", va="bottom", fontsize=5.5)
    ax.axvline(4, color=MUTED, lw=0.4)
    ax.text(4.15, 36.5, r"$t_s\approx 4/|\sigma|$", ha="left", va="bottom", fontsize=5.5)
    ax.text(7.4, 28.5, r"$p=-1\pm1.59j$", ha="right", va="bottom", fontsize=5.5)
    ax.set_title(r"$\theta^*=\pi/4$ (underdamped)")
    ax.set_ylim(27, 61)
    ax.set_yticks([45])
    ax.set_yticklabels([r"$45^\circ$"])

    # θ* = 80°: real poles p1, p2 (overdamped), same initial offset
    a0 = 5 * np.cos(np.radians(80))
    p1, p2 = -1 + np.sqrt(1 - a0), -1 - np.sqrt(1 - a0)
    y = y0 * (p2 * np.exp(p1 * t) - p1 * np.exp(p2 * t)) / (p2 - p1)
    ax = axs[1]
    ax.axhline(80, color=MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.plot(t, 80 + np.degrees(y), color=BLUE)
    ax.axvline(4 / abs(p1), color=MUTED, lw=0.4)
    ax.text(4 / abs(p1) - 0.15, 67.5, r"$t_s\approx 4/0.64$", ha="right", va="bottom", fontsize=5.5)
    ax.text(7.4, 62, r"$p=-0.64,\,-1.36$", ha="right", va="bottom", fontsize=5.5)
    ax.set_title(r"$\theta^*=80^\circ$ (overdamped)")
    ax.set_ylim(61, 86)
    ax.set_yticks([80])
    ax.set_yticklabels([r"$80^\circ$"])

    for ax in axs:
        ax.set_xlim(0, 7.5)
        ax.set_xticks([0, 2, 4, 6])
        ax.set_xlabel(r"$t$", labelpad=0.5)
    fig.subplots_adjust(wspace=0.32)
    save(fig, "pendulum-cases")


def lyap_levelsets():
    """A = [0 1; -2 -0.5] is stable, but circles (P = I) are crossed outward; ALE ellipses never are."""
    k, b = 2.0, 0.5
    A = np.array([[0.0, 1.0], [-k, -b]])
    P_ale = np.array([[(1 + k) / (2 * b) + b / (2 * k), 1 / (2 * k)], [1 / (2 * k), (1 + k) / (2 * b * k)]])
    assert np.allclose(A.T @ P_ale + P_ale @ A, -np.eye(2))

    t = np.linspace(0, 14, 2000)
    x = flow(A, [1.0, 0.0], t)
    pt = flow(A, [1.0, 0.0], [0.3])[:, 0]
    lim = 1.5
    g = np.linspace(-lim, lim, 400)
    X, Y = np.meshgrid(g, g)

    fig, axs = plt.subplots(1, 2, figsize=(2.5, 1.35))
    for ax, P, title in [
        (axs[0], np.eye(2), r"$P=I$: $\dot V>0$ at the dot"),
        (axs[1], P_ale, r"$A^{\mathsf{T}}P+PA=-I$: $\dot V<0$"),
    ]:
        V = P[0, 0] * X**2 + 2 * P[0, 1] * X * Y + P[1, 1] * Y**2
        c = pt @ P @ pt
        ax.contour(X, Y, V, levels=c * np.array([0.25, 0.55, 1.6, 2.4]), colors=AXIS, linewidths=0.4)
        ax.contour(X, Y, V, levels=[c], colors=MUTED, linewidths=0.6)
        ax.plot(x[0], x[1], color=BLUE, lw=0.8)
        normal, velocity = 2 * P @ pt, A @ pt
        n_hat, v_hat = normal / np.linalg.norm(normal), velocity / np.linalg.norm(velocity)
        arrow(ax, pt, 0.4 * n_hat, ORANGE)
        arrow(ax, pt, 0.4 * v_hat, AQUA)
        ax.plot(*pt, marker="o", ms=2.2, color=INK, zorder=4)
        ax.text(*(pt + 0.45 * n_hat + [0, -0.06]), r"$\partial V/\partial \mathbf{x}$", fontsize=5.5, ha="center", va="top")
        ax.text(*(pt + 0.45 * v_hat + [-0.04, 0]), r"$\dot{\mathbf{x}}$", fontsize=5.5, ha="right", va="top")
        ax.set_title(title, fontsize=5.8)
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        bare(ax)
        ax.axhline(0, color=AXIS, lw=0.4, zorder=0)
        ax.axvline(0, color=AXIS, lw=0.4, zorder=0)
    fig.subplots_adjust(wspace=0.08)
    save(fig, "lyap-levelsets")


# lecture-slide crops: name -> (deck, page, x, y, width, height) in 300 dpi pixels
CROPS = {
    "slide-modal-coords": ("lec10", 3, 2700, 640, 2080, 980),
    "slide-uas-ues": ("lec14", 6, 80, 40, 4960, 2360),
    "slide-roa": ("lec17", 1, 2620, 80, 2400, 1832),
}


def slide_crop(name, deck, page, x, y, width, height):
    out = ASSETS / name
    subprocess.run(
        ["pdftoppm", "-r", "300", "-f", str(page), "-l", str(page), "-x", str(x), "-y", str(y),
         "-W", str(width), "-H", str(height), "-png", "-singlefile", str(LECTURES / f"{deck}.pdf"), str(out)],
        check=True,
    )
    # wash out the notebook ruling (light, unsaturated grey) so only the ink and highlights remain
    img = np.asarray(Image.open(f"{out}.png").convert("RGB")).astype(int)
    lo, hi = img.min(axis=2), img.max(axis=2)
    img[(lo > 170) & (hi - lo < 12)] = 255
    Image.fromarray(img.astype(np.uint8)).save(f"{out}.png", optimize=True)


if __name__ == "__main__":
    ASSETS.mkdir(exist_ok=True)
    splane_modes()
    trdet()
    pendulum_cases()
    lyap_levelsets()
    for name, crop in CROPS.items():
        slide_crop(name, *crop)
    print(f"wrote {len(list(ASSETS.iterdir()))} files to {ASSETS}")
