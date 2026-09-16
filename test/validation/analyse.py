#!/usr/bin/env python3
"""Compare the Marlin MuonCVXDDigitiser with its Gaudi port.

Reads the CSV files written by dump_marlin.py, dump_gaudi.py and export_surfaces.py in the work directory,
writes summary.txt (all sub-detectors) and validation_<subdet>_<layer>.png (detailed plots).

With --check, exit with an error if
- the port and Marlin differ in the deterministic part of the chain (cell IDs and times with smearing off), or
- the port with all fixes reverted is not statistically compatible with Marlin.

Usage: analyse.py <work directory> [--subdet VXDBarrel] [--layer 0] [--check]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import IMPLEMENTATIONS, SUBDETECTORS, layer_of  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("work")
parser.add_argument("--subdet", default="VXDBarrel", choices=[s[0] for s in SUBDETECTORS],
                    help="Sub-detector for the detailed plots")
parser.add_argument("--layer", type=int, default=0, help="Layer ID for the detailed plots, -1 for all layers")
parser.add_argument("--check", action="store_true", help="Exit with an error if the validation fails")
args = parser.parse_args()

W = args.work
CONFIGS = ("default", "nosmear")
STYLES = {"marlin": ("#222222", "-"), "port": ("#1f77b4", "-"),
          "revert_cutondeltarays": ("#2ca02c", ":"), "revert_all": ("#ff7f0e", "--")}
KS_MIN_PVALUE = 1e-3

# ---------------------------------------------------------------------------------------------- inputs
surfaces = pd.read_csv(f"{W}/surfaces.csv")
sim = pd.read_csv(f"{W}/sim_hits.csv")
sim["layer"] = layer_of(sim.cellID.values)


def parse_pixels(text):
    if not isinstance(text, str) or not text:
        return []
    return [tuple(float(v) for v in p.split(":")) for p in text.split(";")]


def load(impl, config):
    df = pd.read_csv(f"{W}/hits_{impl}_{config}.csv")
    df = df.merge(surfaces, on=["subdet", "cellID"], how="left")
    if df.ux.isna().any():
        sys.exit(f"hits_{impl}_{config}.csv: cell IDs without surface")
    df["layer"] = layer_of(df.cellID.values)
    d = df[["x", "y", "z"]].values - df[["sim_x", "sim_y", "sim_z"]].values
    p = df[["px", "py", "pz"]].values
    U, V, N = df[["ux", "uy", "uz"]].values, df[["vx", "vy", "vz"]].values, df[["nx", "ny", "nz"]].values
    df["res_u"] = 1e3 * np.einsum("ij,ij->i", d, U)  # um
    df["res_v"] = 1e3 * np.einsum("ij,ij->i", d, V)
    pn = np.abs(np.einsum("ij,ij->i", p, N))
    with np.errstate(invalid="ignore", divide="ignore"):
        df["alpha_u"] = np.degrees(np.arctan(np.einsum("ij,ij->i", p, U) / pn))
        df["alpha_v"] = np.degrees(np.arctan(np.einsum("ij,ij->i", p, V) / pn))
    pixels = df.pixels.map(parse_pixels)
    df["pixset"] = pixels.map(lambda ps: frozenset((x, y) for x, y, _ in ps))
    df["ext_u"] = pixels.map(lambda ps: int(round(max(x for x, _, _ in ps) - min(x for x, _, _ in ps))) + 1 if ps else 0)
    df["ext_v"] = pixels.map(lambda ps: int(round(max(y for _, y, _ in ps) - min(y for _, y, _ in ps))) + 1 if ps else 0)
    df["charge_keV"] = df.edep * 1e6
    df["dt_ps"] = 1e3 * (df.t - df.sim_t)
    df["key"] = list(zip(df.subdet, df.event, df.cellID, df.sim_x.round(5), df.sim_y.round(5), df.sim_z.round(5)))
    return df


impls = [i for i in IMPLEMENTATIONS if all(os.path.exists(f"{W}/hits_{i}_{c}.csv") for c in CONFIGS)]
if "marlin" not in impls or "port" not in impls:
    sys.exit("The Marlin and port outputs are required")
data = {(i, c): load(i, c) for i in impls for c in CONFIGS}


def select(df, subdet, layer=-1):
    df = df[df.subdet == subdet]
    return df if layer < 0 else df[df.layer == layer]


# ---------------------------------------------------------------------------------------------- summary
core = lambda v: (np.percentile(v, 75) - np.percentile(v, 25)) / 1.349  # noqa: E731
failures = []


def compare_hits(ref, other):
    both = ref.merge(other, on="key", suffixes=("_r", "_o"))
    dpos = 1e3 * np.sqrt((both.x_r - both.x_o) ** 2 + (both.y_r - both.y_o) ** 2 + (both.z_r - both.z_o) ** 2)
    overlap = np.array([len(a & b) / len(a | b) if (a | b) else 1. for a, b in zip(both.pixset_r, both.pixset_o)])
    return both, dpos, overlap


def summary_block(title, subdet, layer):
    nsim = len(select(sim, subdet, layer))
    if nsim == 0:
        return [f"--- {title}: no SimTrackerHits"]
    with_momentum = select(sim, subdet, layer).has_momentum.mean()
    lines = [f"--- {title}: {nsim} SimTrackerHits ({with_momentum:.1%} with momentum and path length)"]
    for config in CONFIGS:
        ref = select(data[("marlin", config)], subdet, layer)
        tail_charge = 2.5 * ref.charge_keV.median()
        lines.append(f"  [{config}] {'':36s} {'eff':>6s}  {'res u core/tail':>21s}  {'res v core/tail':>21s}  "
                     f"{'size':>12s}  {f'charge med/>{tail_charge:.0f} keV':>24s}  {'t core':>7s}")
        for impl in impls:
            df = select(data[(impl, config)], subdet, layer)
            ks = {q: stats.ks_2samp(ref[q], df[q]).pvalue for q in ("res_u", "res_v", "npix", "charge_keV")}
            tag = (lambda q: "") if impl == "marlin" else (lambda q: f" ({ks[q]:.2f})")
            res = lambda q: f"{core(df[q]):5.2f}/{np.mean(np.abs(df[q]) > 15):.3f}{tag(q)}"  # noqa: E731
            charge = f"{df.charge_keV.median():5.1f}/{np.mean(df.charge_keV > tail_charge):.3f}{tag('charge_keV')}"
            lines.append(f"  {IMPLEMENTATIONS[impl]:45s} {len(df) / nsim:6.4f}  {res('res_u'):>21s}  "
                         f"{res('res_v'):>21s}  {df.npix.mean():5.2f}{tag('npix'):>7s}  {charge:>24s}  "
                         f"{core(df.dt_ps):7.1f}")
            if impl == "revert_all" and layer < 0:
                eff_err = np.sqrt(len(ref) * (1 - len(ref) / nsim)) / nsim
                if abs(len(df) - len(ref)) / nsim > max(5 * eff_err, 1e-3):
                    failures.append(f"{subdet} [{config}]: efficiency of the port with all fixes reverted differs from Marlin")
                for q, p in ks.items():
                    if p < KS_MIN_PVALUE:
                        failures.append(f"{subdet} [{config}]: {q} of the port with all fixes reverted differs from Marlin (KS p={p:.1e})")
    lines.append("  [nosmear, hit by hit]")
    for impl in impls[1:]:
        both, dpos, overlap = compare_hits(select(data[("marlin", "nosmear")], subdet, layer),
                                           select(data[(impl, "nosmear")], subdet, layer))
        same_cell, same_time = np.mean(both.cellID_r == both.cellID_o), np.mean(np.isclose(both.t_r, both.t_o))
        lines.append(f"  {IMPLEMENTATIONS[impl]:45s} matched {len(both):6d}  same cellID {same_cell:.4f}  "
                     f"same time {same_time:.4f}  |dpos| < 1 um {np.mean(dpos < 1):.3f}  "
                     f"identical pixels {np.mean(overlap == 1):.3f}  pixel overlap {np.mean(overlap):.3f}")
        if layer < 0 and (same_cell < 1 or same_time < 1):
            failures.append(f"{subdet}: {IMPLEMENTATIONS[impl]} and Marlin differ in cell IDs or times without smearing")
    return lines


lines = ["Validation of the Gaudi port of MuonCVXDDigitiser against the Marlin processor",
         "residuals: core width IQR/1.349 [um] / fraction beyond 15 um; t core [ps]; "
         "KS p-values against Marlin in brackets", ""]
for prefix, *_ in SUBDETECTORS:
    lines += summary_block(prefix, prefix, -1) + [""]
layer_title = f"{args.subdet} layer {args.layer}" if args.layer >= 0 else f"{args.subdet} all layers"
if args.layer >= 0:
    lines += summary_block(layer_title, args.subdet, args.layer) + [""]
lines += ["Checks: " + ("all passed" if not failures else f"{len(failures)} FAILED")] + [f"  FAILED: {f}" for f in failures]
report = "\n".join(lines)
print(report)
open(f"{W}/summary.txt", "w").write(report + "\n")

# ---------------------------------------------------------------------------------------------- plots
plt.rcParams.update({"font.size": 9})
fig, axes = plt.subplots(3, 4, figsize=(16, 11))
plot_data = {k: select(v, args.subdet, args.layer) for k, v in data.items()}


def hist(ax, q, bins, xlabel, logy=False):
    for impl in impls:
        color, ls = STYLES[impl]
        v = plot_data[(impl, "default")][q]
        ax.hist(np.clip(v, bins[0], bins[-1]), bins=bins, histtype="step", color=color, ls=ls, lw=1.4,
                density=True, label=IMPLEMENTATIONS[impl])
    ax.set_xlabel(xlabel)
    if logy:
        ax.set_yscale("log")


def profile(ax, x, y, bins, xlabel, ylabel, fn=np.mean):
    centers = 0.5 * (bins[1:] + bins[:-1])
    for impl in impls:
        color, ls = STYLES[impl]
        df = plot_data[(impl, "default")]
        idx = np.digitize(df[x], bins) - 1
        values = [fn(df[y][idx == i]) if np.sum(idx == i) > 20 else np.nan for i in range(len(centers))]
        ax.plot(centers, values, color=color, ls=ls, marker="o", ms=3, label=IMPLEMENTATIONS[impl])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


rms = lambda v: np.sqrt(np.mean(np.square(v)))  # noqa: E731
qmax = 5 * plot_data[("marlin", "default")].charge_keV.median()
hist(axes[0, 0], "res_u", np.linspace(-40, 40, 81), "residual u [um]", logy=True)
hist(axes[0, 1], "res_v", np.linspace(-40, 40, 81), "residual v [um]", logy=True)
hist(axes[0, 2], "npix", np.arange(0.5, 20.5, 1), "cluster size [pixels]")
hist(axes[0, 3], "charge_keV", np.linspace(0, qmax, 51), "cluster charge [keV]", logy=True)
hist(axes[1, 0], "ext_u", np.arange(0.5, 10.5, 1), "cluster extent u [pixels]")
hist(axes[1, 1], "ext_v", np.arange(0.5, 16.5, 1), "cluster extent v [pixels]")
hist(axes[1, 2], "dt_ps", np.linspace(-150, 150, 61), "t - t_sim [ps]")
profile(axes[1, 3], "alpha_v", "npix", np.linspace(-60, 60, 13), "incidence angle along v [deg]",
        "mean cluster size [pixels]")
profile(axes[2, 0], "alpha_u", "res_u", np.linspace(-40, 40, 9), "incidence angle along u [deg]",
        "RMS residual u [um]", fn=rms)
profile(axes[2, 1], "alpha_v", "res_v", np.linspace(-60, 60, 13), "incidence angle along v [deg]",
        "RMS residual v [um]", fn=rms)
for impl in impls[1:]:
    color, ls = STYLES[impl]
    _, dpos, overlap = compare_hits(plot_data[("marlin", "nosmear")], plot_data[(impl, "nosmear")])
    axes[2, 2].hist(np.clip(dpos, 0, 30), bins=np.linspace(0, 30, 61), histtype="step", color=color, ls=ls, lw=1.4,
                    density=True, label=IMPLEMENTATIONS[impl])
    axes[2, 3].hist(overlap, bins=np.linspace(0, 1, 21), histtype="step", color=color, ls=ls, lw=1.4, density=True,
                    label=IMPLEMENTATIONS[impl])
axes[2, 2].set_yscale("log")
axes[2, 2].set_xlabel("|position - Marlin position| [um], no smearing")
axes[2, 3].set_xlabel("fired pixel overlap with Marlin (Jaccard), no smearing")
for ax in axes.flat:
    ax.grid(alpha=0.3)
axes[0, 0].legend(loc="upper left", fontsize=8)
axes[2, 2].legend(loc="upper right", fontsize=8)
nsim_plot = len(select(sim, args.subdet, args.layer))
fig.suptitle(f"MuonCVXDDigitiser, Marlin vs Gaudi port: {layer_title} ({nsim_plot} SimTrackerHits); "
             "default settings unless stated", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
suffix = f"layer{args.layer}" if args.layer >= 0 else "all"
fig.savefig(f"{W}/validation_{args.subdet}_{suffix}.png", dpi=130)
print(f"Plots written to {W}/validation_{args.subdet}_{suffix}.png")

if args.check and failures:
    sys.exit(1)
