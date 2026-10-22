#!/usr/bin/env python3
"""Cluster shapes of MuonCVXDDigitiser from the fired pixels (StoreFiredPixels), with square and long tracker
pixels, and resolutions compared with DDPlanarDigi.

Reads hits_<name>.csv, sim_hits.csv and surfaces.csv from the work directory, writes clusters_summary.txt,
clusters_shapes.png, clusters_resolution_v.png and clusters_examples.png.

Usage: clusters.py <work directory>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

W = sys.argv[1]
SUBDETS = ["VXDBarrel", "VXDEndcap", "ITBarrel", "ITEndcap", "OTBarrel", "OTEndcap"]
# MAIA_v0 sensor thickness (mm) and MuonCVXDDigitiser TanLorentz
THICKNESS = {"VXDBarrel": 0.05, "VXDEndcap": 0.05, "ITBarrel": 0.1, "ITEndcap": 0.1, "OTBarrel": 0.1, "OTEndcap": 0.1}
TAN_LORENTZ = 0.8
PITCH_U = 0.025
# name: (label, colour, line style, tracker pixel size along v in mm)
IMPLS = {"ddplanar_fixed": ("DDPlanarDigi (parametric)", "#d62728", "-", None),
         "cvxd_square": ("MuonCVXDDigitiser, 25x25 um pixels", "#1f77b4", "-", 0.025),
         "cvxd_long": ("MuonCVXDDigitiser, trackers 25x312 um pixels", "#2ca02c", "--", 0.312)}
REALISTIC = ["cvxd_square", "cvxd_long"]

surf = pd.read_csv(f"{W}/surfaces.csv")
sim = pd.read_csv(f"{W}/sim_hits.csv")


def parse(text):
    if not isinstance(text, str) or not text:
        return []
    return [tuple(float(x) for x in p.split(":")) for p in text.split(";")]


def pitch_v(impl, sd):
    return 0.025 if sd.startswith("VXD") else IMPLS[impl][3]


def load(impl):
    df = pd.read_csv(f"{W}/hits_{impl}.csv").merge(surf, on=["subdet", "cellID"], how="left")
    d = df[["x", "y", "z"]].values - df[["sim_x", "sim_y", "sim_z"]].values
    U, V, N = df[["ux", "uy", "uz"]].values, df[["vx", "vy", "vz"]].values, df[["nx", "ny", "nz"]].values
    df["res_u"] = 1e3 * np.einsum("ij,ij->i", d, U)
    df["res_v"] = 1e3 * np.einsum("ij,ij->i", d, V)
    df["pull_u"] = df.res_u / (1e3 * df.du)
    df["pull_v"] = df.res_v / (1e3 * df.dv)
    p = df[["px", "py", "pz"]].values
    # signed: the Lorentz drift shears the charge in a fixed local direction, whichever side the particle enters
    pn = np.einsum("ij,ij->i", p, N)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["tan_u"] = np.einsum("ij,ij->i", p, U) / pn
        df["tan_v"] = np.einsum("ij,ij->i", p, V) / pn
    df["alpha_u"] = np.degrees(np.arctan(df.tan_u))
    df["alpha_v"] = np.degrees(np.arctan(df.tan_v))
    if df.pixels.notna().any():
        pix = df.pixels.map(parse)
        df["npix"] = pix.map(len)
        df["ext_u"] = pix.map(lambda ps: int(round(max(p[0] for p in ps) - min(p[0] for p in ps))) + 1 if ps else 0)
        df["ext_v"] = pix.map(lambda ps: int(round(max(p[1] for p in ps) - min(p[1] for p in ps))) + 1 if ps else 0)
        df["dt_pix_ps"] = pix.map(lambda ps: 1e3 * (max(p[3] for p in ps) - min(p[3] for p in ps)) if ps else np.nan)
        df["pix_list"] = pix
    df["key"] = list(zip(df.subdet, df.event, df.cellID, df.sim_x.round(5), df.sim_y.round(5), df.sim_z.round(5)))
    return df


data = {impl: load(impl) for impl in IMPLS}
for impl in REALISTIC:
    if "npix" not in data[impl]:
        sys.exit(f"hits_{impl}.csv has no fired pixels: run the realistic digitisation with --store-pixels")
core = lambda v: (np.percentile(v, 75) - np.percentile(v, 25)) / 1.349 if len(v) else np.nan  # noqa: E731

# ------------------------------------------------------------------------------------------- summary
lines = ["Cluster shapes and resolutions, MAIAConfig tracker digitisation, "
         f"{len(sim)} SimTrackerHits in {sim.event.nunique()} events", "",
         f"{'':10s} {'':16s} {'eff':>6s} {'size':>5s} {'ext u':>6s} {'ext v':>6s} {'res u':>6s} {'res v':>7s} "
         f"{'du':>6s} {'dv':>6s} {'pull u':>6s} {'pull v':>6s} {'pix t spread':>12s}"]
for sd in SUBDETS:
    nsim = (sim.subdet == sd).sum()
    for impl in IMPLS:
        df = data[impl][data[impl].subdet == sd]
        has_pix = "npix" in df
        lines.append(
            f"{sd:10s} {impl:16s} {df.key.nunique() / nsim:6.3f} "
            + (f"{df.npix.mean():5.2f} {df.ext_u.mean():6.2f} {df.ext_v.mean():6.2f} " if has_pix
               else f"{'-':>5s} {'-':>6s} {'-':>6s} ")
            + f"{core(df.res_u):6.2f} {core(df.res_v):7.2f} {1e3 * df.du.median():6.2f} {1e3 * df.dv.median():6.2f} "
              f"{core(df.pull_u):6.2f} {core(df.pull_v):6.2f} "
            + (f"{df.dt_pix_ps.mean():9.1f} ps" if has_pix else f"{'-':>12s}"))
    lines.append("")
lines.append("res: core width IQR/1.349 [um]; size: fired pixels per cluster; ext: cluster extent in pixels; "
             "pix t spread: mean (max - min) of the pixel times in a cluster")
report = "\n".join(lines)
print(report)
open(f"{W}/clusters_summary.txt", "w").write(report + "\n")

# ------------------------------------------------------------------------------------------- figure 1: shapes
plt.rcParams.update({"font.size": 8})
fig, axes = plt.subplots(5, 6, figsize=(20, 15))
for j, sd in enumerate(SUBDETS):
    for impl in REALISTIC:
        if impl == "cvxd_long" and sd.startswith("VXD"):
            continue  # same configuration as cvxd_square for the vertex detector
        label, color, ls, _ = IMPLS[impl]
        df = data[impl][data[impl].subdet == sd]
        kw = dict(histtype="step", color=color, ls=ls, lw=1.3, density=True, label=label)
        axes[0, j].hist(df.npix.clip(upper=20), bins=np.arange(0.5, 21.5, 1), **kw)
        axes[1, j].hist(df.ext_u.clip(upper=10), bins=np.arange(0.5, 11.5, 1), **kw)
        axes[2, j].hist(df.ext_v.clip(upper=20), bins=np.arange(0.5, 21.5, 1), **kw)
        # extent vs incidence angle, with the geometric expectation
        for row, (ang, ext, pitch, lorentz) in {3: ("tan_u", "ext_u", PITCH_U, TAN_LORENTZ),
                                                4: ("tan_v", "ext_v", pitch_v(impl, sd), 0.)}.items():
            bins = np.linspace(-2.5, 2.5, 21) if row == 3 else np.linspace(-4, 4, 21)
            idx = np.digitize(df[ang], bins) - 1
            centers = 0.5 * (bins[1:] + bins[:-1])
            mean = [df[ext][idx == i].mean() if np.sum(idx == i) > 15 else np.nan for i in range(len(centers))]
            axes[row, j].plot(centers, mean, color=color, ls=ls, marker="o", ms=3, label=label)
            expect = THICKNESS[sd] * np.abs(centers - lorentz) / pitch + 1
            axes[row, j].plot(centers, expect, color=color, ls=":", lw=1,
                              label="thickness x |tan - tan_L| / pitch + 1" if row == 3 else "thickness x |tan| / pitch + 1")
    axes[0, j].set_title(sd)
    axes[0, j].set_xlabel("cluster size [pixels]")
    axes[1, j].set_xlabel("cluster extent u [pixels]")
    axes[2, j].set_xlabel("cluster extent v [pixels]")
    axes[3, j].set_xlabel("tan(incidence angle along u)")
    axes[3, j].set_ylabel("mean extent u [pixels]")
    axes[4, j].set_xlabel("tan(incidence angle along v)")
    axes[4, j].set_ylabel("mean extent v [pixels]")
for ax in axes.flat:
    ax.grid(alpha=0.3)
axes[0, 2].legend(fontsize=7, loc="upper right")
axes[3, 2].legend(fontsize=6, loc="upper center")
fig.suptitle("MuonCVXDDigitiser cluster shapes from the fired pixels: square 25x25 um pixels and long 25x312 um tracker "
             f"pixels (vertex unchanged); expectation includes the Lorentz drift (tan = {TAN_LORENTZ}) along u",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{W}/clusters_shapes.png", dpi=110)

# ------------------------------------------------------------------------------------------- figure 2: resolution
fig, axes = plt.subplots(3, 6, figsize=(20, 9))
for j, sd in enumerate(SUBDETS):
    for impl, (label, color, ls, _) in IMPLS.items():
        if impl == "cvxd_long" and sd.startswith("VXD"):
            continue
        df = data[impl][data[impl].subdet == sd]
        kw = dict(histtype="step", color=color, ls=ls, lw=1.3, density=True, label=label)
        lim = 40 if sd.startswith("VXD") else 400
        axes[0, j].hist(df.res_v.clip(-lim, lim), bins=np.linspace(-lim, lim, 81), **kw)
        axes[1, j].hist(df.pull_v.clip(-6, 6), bins=np.linspace(-6, 6, 61), **kw)
        if "npix" in df and len(df):
            q = np.concatenate([[p[2] for p in ps] for ps in df.pix_list])
            axes[2, j].hist(np.clip(q, 0, 4e4), bins=np.linspace(0, 4e4, 81), **kw)
    axes[0, j].set_title(sd)
    axes[0, j].set_xlabel("residual v [um]")
    axes[0, j].set_yscale("log")
    axes[1, j].set_xlabel("pull v = residual v / dv")
    axes[1, j].set_yscale("log")
    x = np.linspace(-6, 6, 200)
    axes[1, j].plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color="gray", ls=":", lw=0.8)
    axes[2, j].set_xlabel("charge per fired pixel [electrons]")
    axes[2, j].set_yscale("log")
for ax in axes.flat:
    ax.grid(alpha=0.3)
axes[0, 2].legend(fontsize=7, loc="upper left")
fig.suptitle("Resolution along v and pixel charge: DDPlanarDigi (7 x 90 um in the trackers) vs MuonCVXDDigitiser with "
             "square and long tracker pixels", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{W}/clusters_resolution_v.png", dpi=110)

# ------------------------------------------------------------------------------------------- figure 3: examples
examples = [("VXDBarrel", "cvxd_square"), ("ITBarrel", "cvxd_square"), ("ITBarrel", "cvxd_long"),
            ("OTEndcap", "cvxd_long")]
angle_bins = [(0, 10), (20, 35), (45, 60), (65, 80)]
fig, axes = plt.subplots(len(examples), len(angle_bins), figsize=(16, 16))
rng = np.random.default_rng(1)
for i, (sd, impl) in enumerate(examples):
    df = data[impl][data[impl].subdet == sd]
    for k, (lo, hi) in enumerate(angle_bins):
        ax = axes[i, k]
        sel = df[(df.alpha_v.abs() >= lo) & (df.alpha_v.abs() < hi) & (df.npix > 0)]
        if sel.empty:
            ax.set_axis_off()
            continue
        row = sel.iloc[rng.integers(len(sel))]
        ps = row.pix_list
        us, vs, qs = np.array([p[0] for p in ps]), np.array([p[1] for p in ps]), np.array([p[2] for p in ps])
        u0, v0 = us.min(), vs.min()
        nu, nv = int(round(us.max() - u0)) + 1, int(round(vs.max() - v0)) + 1
        grid = np.full((nv, nu), np.nan)
        for uu, vv, qq in zip(us, vs, qs):
            grid[int(round(vv - v0)), int(round(uu - u0))] = qq
        # physical aspect ratio of the pixels
        ax.imshow(grid, origin="lower", cmap="viridis", aspect=pitch_v(impl, sd) / PITCH_U,
                  extent=(-0.5, nu - 0.5, -0.5, nv - 0.5))
        for (vv, uu), qq in np.ndenumerate(grid):
            if not np.isnan(qq):
                ax.text(uu, vv, f"{qq / 1e3:.1f}k", ha="center", va="center", fontsize=6, color="white")
        ax.set_title(f"{sd} {'long' if impl == 'cvxd_long' else 'square'} pixels, |alpha_v| = {abs(row.alpha_v):.0f} deg, "
                     f"alpha_u = {row.alpha_u:.0f} deg\nres u = {row.res_u:.1f} um, res v = {row.res_v:.1f} um",
                     fontsize=7)
        ax.set_xlabel("u [pixel]")
        ax.set_ylabel("v [pixel]")
fig.suptitle("Example clusters (charge per pixel in electrons; pixels drawn with their physical aspect ratio)",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{W}/clusters_examples.png", dpi=110)
print(f"Plots written to {W}/clusters_shapes.png, clusters_resolution_v.png and clusters_examples.png")
