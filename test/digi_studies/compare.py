#!/usr/bin/env python3
"""Compare the tracker hits of DDPlanarDigi (parametric) and MuonCVXDDigitiser (realistic) as configured in
MAIAConfig: residuals and pulls in the local frame of the sensors, time convention, energy, efficiency and
the EDM fields each digitiser fills.

Reads hits_<name>.csv, sim_hits.csv and surfaces.csv from the work directory, writes compare_summary.txt,
compare_resolution.png and compare_time_energy.png.

Usage: compare.py <work directory>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

W = sys.argv[1]
C_MM_NS = 299.792458
SUBDETS = ["VXDBarrel", "VXDEndcap", "ITBarrel", "ITEndcap", "OTBarrel", "OTEndcap"]
ALL_IMPLS = {"ddplanar": ("DDPlanarDigi, MAIAConfig as is", "#bbbbbb", "--"),
             "ddplanar_fixed": ("DDPlanarDigi, time window corrected for r/c", "#d62728", "-"),
             "cvxd_long": ("MuonCVXDDigitiser, MAIAConfig as is", "#1f77b4", "-")}
# The distributions compare the parametric configuration, with the time window corrected, with the realistic one
IMPLS = {k: ALL_IMPLS[k] for k in ("ddplanar_fixed", "cvxd_long")}

surf = pd.read_csv(f"{W}/surfaces.csv")
sim = pd.read_csv(f"{W}/sim_hits.csv")
sim["key"] = list(zip(sim.subdet, sim.event, sim.cellID, sim.x.round(5), sim.y.round(5), sim.z.round(5)))
sim["t_corr"] = sim.t - np.sqrt(sim.x ** 2 + sim.y ** 2 + sim.z ** 2) / C_MM_NS


def load(impl):
    df = pd.read_csv(f"{W}/hits_{impl}.csv").merge(surf, on=["subdet", "cellID"], how="left")
    if df.ux.isna().any():
        sys.exit(f"hits_{impl}.csv: cell IDs without surface")
    d = df[["x", "y", "z"]].values - df[["sim_x", "sim_y", "sim_z"]].values
    df["res_u"] = 1e3 * np.einsum("ij,ij->i", d, df[["ux", "uy", "uz"]].values)  # um
    df["res_v"] = 1e3 * np.einsum("ij,ij->i", d, df[["vx", "vy", "vz"]].values)
    df["pull_u"] = df.res_u / (1e3 * df.du)
    df["pull_v"] = df.res_v / (1e3 * df.dv)
    r_sim = np.sqrt(df.sim_x ** 2 + df.sim_y ** 2 + df.sim_z ** 2)
    df["dt_ps"] = 1e3 * (df.t - df.sim_t)
    df["dt_corr_ps"] = 1e3 * (df.t - (df.sim_t - r_sim / C_MM_NS))
    df["edep_ratio"] = df.edep / df.sim_edep
    df["key"] = list(zip(df.subdet, df.event, df.cellID, df.sim_x.round(5), df.sim_y.round(5), df.sim_z.round(5)))
    df["has_cov"] = df[[f"cov{i}" for i in range(6)]].abs().sum(axis=1) > 0
    return df


data = {impl: load(impl) for impl in ALL_IMPLS}
core = lambda v: (np.percentile(v, 75) - np.percentile(v, 25)) / 1.349 if len(v) else np.nan  # noqa: E731

# ------------------------------------------------------------------------------------------- summary
lines = ["DDPlanarDigi vs MuonCVXDDigitiser as configured in MAIAConfig",
         f"{len(sim)} SimTrackerHits in {sim.event.nunique()} events", "",
         "Residuals: core width IQR/1.349 [um] / RMS [um]; pulls: core width; times: median / core [ps]", "",
         f"{'':10s} {'':14s} {'eff':>6s} {'res u':>13s} {'res v':>13s} {'du':>6s} {'dv':>6s} {'pull u':>6s} "
         f"{'pull v':>6s} {'t-t_sim':>13s} {'t-(t_sim-r/c)':>14s} {'E/E_sim':>13s} {'cov':>4s}"]
for sd in SUBDETS:
    nsim = (sim.subdet == sd).sum()
    for impl in ALL_IMPLS:
        df = data[impl][data[impl].subdet == sd]
        if df.empty:
            lines.append(f"{sd:10s} {impl:14s} {0:6.3f}  no hits")
            continue
        lines.append(
            f"{sd:10s} {impl:14s} {df.key.nunique() / nsim:6.3f} "
            f"{core(df.res_u):5.2f}/{df.res_u.std():6.2f} {core(df.res_v):5.2f}/{df.res_v.std():6.2f} "
            f"{1e3 * df.du.median():6.2f} {1e3 * df.dv.median():6.2f} {core(df.pull_u):6.2f} {core(df.pull_v):6.2f} "
            f"{df.dt_ps.median():6.1f}/{core(df.dt_ps):6.1f} {df.dt_corr_ps.median():7.1f}/{core(df.dt_corr_ps):6.1f} "
            f"{df.edep_ratio.median():5.3f}/{core(df.edep_ratio):6.3f} {df.has_cov.mean():4.2f}")
    lines.append("")

lines.append("EDM fields of the TrackerHitPlanes (all sub-detectors)")
for impl, (label, _, _) in ALL_IMPLS.items():
    df = data[impl]
    lines.append(f"  {label}:")
    lines.append(f"    du values [um]: {sorted(set(np.round(1e3 * df.du, 2)))}   "
                 f"dv values [um]: {sorted(set(np.round(1e3 * df.dv, 2)))}")
    lines.append(f"    covMatrix filled: {df.has_cov.mean():.3f}   type: {sorted(set(df.type))}   "
                 f"quality: {sorted(set(df.quality))}   edepError: {sorted(set(df.edepError))[:3]}")
    lines.append(f"    cellID equal to the SimTrackerHit: {np.mean(df.cellID == df.reco_cellID):.3f}   "
                 f"links per SimTrackerHit max: {df.nlinks_sim.max()}")
    lines.append("    distinct times per sub-detector (time discretisation): "
                 + ", ".join(f"{sd} {df[df.subdet == sd].t.round(6).nunique()}" for sd in SUBDETS))
report = "\n".join(lines)
print(report)
open(f"{W}/compare_summary.txt", "w").write(report + "\n")

# ------------------------------------------------------------------------------------------- plots
plt.rcParams.update({"font.size": 8})


def hist(ax, q, bins, sd, logy=False):
    for impl, (label, color, ls) in IMPLS.items():
        v = data[impl][data[impl].subdet == sd][q]
        ax.hist(np.clip(v, bins[0], bins[-1]), bins=bins, histtype="step", color=color, ls=ls, lw=1.3,
                density=True, label=label)
    if logy:
        ax.set_yscale("log")


# Figure 1: spatial residuals and pulls
fig, axes = plt.subplots(4, 6, figsize=(20, 12))
for j, sd in enumerate(SUBDETS):
    hist(axes[0, j], "res_u", np.linspace(-40, 40, 81), sd, logy=True)
    lim_v = 40 if sd.startswith("VXD") else 300
    hist(axes[1, j], "res_v", np.linspace(-lim_v, lim_v, 81), sd, logy=True)
    hist(axes[2, j], "pull_u", np.linspace(-6, 6, 61), sd, logy=True)
    hist(axes[3, j], "pull_v", np.linspace(-6, 6, 61), sd, logy=True)
    axes[0, j].set_title(sd)
    axes[0, j].set_xlabel("residual u [um]")
    axes[1, j].set_xlabel("residual v [um]")
    axes[2, j].set_xlabel("pull u = residual u / du")
    axes[3, j].set_xlabel("pull v = residual v / dv")
    for i in (2, 3):
        x = np.linspace(-6, 6, 200)
        axes[i, j].plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color="gray", lw=0.8, ls=":",
                        label="unit Gaussian")
for ax in axes.flat:
    ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7, loc="upper left")
axes[2, 0].legend(fontsize=7, loc="upper left")
fig.suptitle("Tracker hits in MAIAConfig: DDPlanarDigi vs MuonCVXDDigitiser, local residuals with respect to the "
             "SimTrackerHit and pulls with the stored uncertainties du, dv (used by ACTS)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{W}/compare_resolution.png", dpi=110)

# Figure 2: time, energy, efficiency
fig, axes = plt.subplots(3, 6, figsize=(20, 9))
for j, sd in enumerate(SUBDETS):
    for impl, (label, color, ls) in IMPLS.items():
        df = data[impl][data[impl].subdet == sd]
        axes[0, j].hist(np.clip(df.dt_ps, -6000, 400), bins=np.linspace(-6000, 400, 129), histtype="step",
                        color=color, lw=1.3, density=True, label=f"{label}: t - t_sim")
        axes[0, j].hist(np.clip(df.dt_corr_ps, -6000, 400), bins=np.linspace(-6000, 400, 129), histtype="step",
                        color=color, lw=1.0, ls=":", density=True, label=f"{label}: t - (t_sim - r/c)")
        axes[1, j].hist(np.clip(df.edep_ratio, 0, 3), bins=np.linspace(0, 3, 61), histtype="step", color=color,
                        lw=1.3, density=True, label=label)
    # efficiency vs propagation-corrected SimTrackerHit time
    s_sd = sim[sim.subdet == sd]
    bins = np.linspace(-0.3, 1.0, 27)
    idx = np.digitize(s_sd.t_corr, bins) - 1
    centers = 0.5 * (bins[1:] + bins[:-1])
    for impl, (label, color, ls) in ALL_IMPLS.items():
        found = s_sd.key.isin(set(data[impl][data[impl].subdet == sd].key))
        eff = [found[idx == i].mean() if np.sum(idx == i) > 10 else np.nan for i in range(len(centers))]
        axes[2, j].plot(centers, eff, color=color, ls=ls, marker="o", ms=3, label=label)
    axes[0, j].set_yscale("log")
    axes[0, j].set_title(sd)
    axes[0, j].set_xlabel("t_hit - t_sim [ps]")
    axes[1, j].set_xlabel("EDep_hit / EDep_sim")
    axes[1, j].set_yscale("log")
    axes[2, j].set_xlabel("SimTrackerHit t - r/c [ns]")
    axes[2, j].set_ylabel("efficiency")
    axes[2, j].set_ylim(0, 1.05)
for ax in axes.flat:
    ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=6, loc="upper left")
axes[2, 0].legend(fontsize=7, loc="lower left")
fig.suptitle("Tracker hits in MAIAConfig: DDPlanarDigi vs MuonCVXDDigitiser, time convention, energy and efficiency",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{W}/compare_time_energy.png", dpi=110)
print(f"Plots written to {W}/compare_resolution.png and {W}/compare_time_energy.png")
