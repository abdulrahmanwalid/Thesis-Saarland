def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
warnings.filterwarnings("ignore")

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                   # headless (cluster) — saves figures, no display
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms

# ═══════════════════════════════════════════════════════════════════════════
# RMSD and MM-GBSA analysis of AD51-docked MD trajectories.
#
# Expected under BASE:
#     complex.prmtop                              solvated topology (matches the DCD)
#     <replica dir>/trajectory.dcd                production run
#     <replica dir>/trajectory_continuation.dcd   optional second leg
#     <replica dir>/trajectory_imaged.dcd         optional, from cpptraj autoimage
#
# Additionally required when RUN_MMGBSA = True:
#     complex_stripped.prmtop   solute only, receptor atoms BEFORE ligand atoms
#     receptor.prmtop
#     ligand.prmtop
#
# Frame 0 is the alignment reference (RMSD = 0 by definition) and is excluded
# from every reported RMSD mean. MM-GBSA keeps frame 0: there it is a real value.
# ═══════════════════════════════════════════════════════════════════════════

# ─────────────────────────── CONFIG — edit this block ──────────────────────
NAME       = "SOX Mg Pose 1"            # title shown on figures
TAG        = "sox_mg_pose1"             # prefix for output file names
BASE       = "/path/to/system"          # folder holding the prmtops
REPS       = [1, 2, 3, 4, 5]            # finished replicas only
REP_DIRS   = {r: f"{BASE}/Mg_Dock_{r}_3POV" for r in REPS}   # edit to your layout
PROT_COLOR = "tab:blue"                 # SOX blue | EBV orange | UL12 green | UL98 red
RUN_MMGBSA = True                       # False for RMSD only

LIG   = "LIG"                           # ligand residue name
DT_NS = 0.01                            # 10 ps per frame
# ───────────────────────────────────────────────────────────────────────────

TOP   = f"{BASE}/complex.prmtop"
CACHE = f"{BASE}/cache";  os.makedirs(CACHE, exist_ok=True)
GRAPH = f"{BASE}/graphs"; os.makedirs(GRAPH, exist_ok=True)

PARTS = {"Backbone": "steelblue", "Protein": PROT_COLOR, "Ligand": "crimson"}

# sized to REPS, so adding a replica can never desync the colour list
REP_PALETTE = ["tab:gray", "tab:purple", "tab:brown", "tab:olive",
               "tab:cyan", "tab:pink", "tab:blue", "tab:orange"]
REP_C = {r: REP_PALETTE[i % len(REP_PALETTE)] for i, r in enumerate(REPS)}

CONT_NAME   = "trajectory_continuation.dcd"
IMAGED_NAME = "trajectory_imaged.dcd"
N_SAMPLE_GPU, N_SAMPLE_CPU = 1500, 500


def rep_input(rep):
    """Trajectory file list for one replica, plus whether it is the imaged copy."""
    d   = REP_DIRS[rep]
    img = f"{d}/{IMAGED_NAME}"
    if os.path.exists(img):
        raw = [f"{d}/trajectory.dcd"]
        if os.path.exists(f"{d}/{CONT_NAME}"):
            raw.append(f"{d}/{CONT_NAME}")
        n_img = mda.Universe(TOP, img).trajectory.n_frames
        n_raw = mda.Universe(TOP, raw).trajectory.n_frames
        if n_img != n_raw:
            raise RuntimeError(
                f"Rep {rep}: imaged has {n_img} frames, raw has {n_raw} — reimage it")
        return [img], True
    legs = [f"{d}/trajectory.dcd"]
    if os.path.exists(f"{d}/{CONT_NAME}"):
        legs.append(f"{d}/{CONT_NAME}")
    return legs, False


inputs = {r: rep_input(r) for r in REPS}
sfx    = lambda imaged: "_imaged" if imaged else ""

# ════════════════════════════════ 1. RMSD ══════════════════════════════════
# backbone / protein / ligand, fitted on protein backbone, reference = frame 0
rmsd = {}
for rep in REPS:
    legs, imaged = inputs[rep]
    f = f"{CACHE}/rmsd_rep{rep}{sfx(imaged)}.npy"
    if os.path.exists(f):
        r = np.load(f)
        print(f"\nRep {rep}: cache{' (imaged)' if imaged else ''}")
    else:
        print(f"\nRep {rep}: {', '.join(os.path.basename(x) for x in legs)}")
        u = mda.Universe(TOP, legs)
        r = rms.RMSD(u, u, select="backbone",
                     groupselections=["protein", f"resname {LIG}"],
                     ref_frame=0).run(verbose=True).results.rmsd[:, 2:]
        np.save(f, r)
    rmsd[rep] = r

    print(f"  {len(r)} frames ({len(r) * DT_NS:.1f} ns)", flush=True)
    for j, part in enumerate(PARTS):
        c = r[1:, j]                                        # frame 0 excluded
        print(f"  {part:<9} min {c.min():6.2f} | max {c.max():6.2f} | "
              f"mean {c.mean():6.2f} +- {c.std():.2f} | "
              f"max frame-to-frame jump {np.abs(np.diff(c)).max():6.2f} A", flush=True)

lengths = {r: len(a) for r, a in rmsd.items()}
if max(lengths.values()) - min(lengths.values()) > 0.02 * max(lengths.values()):
    print("\n*** WARNING: replica lengths differ by more than 2% ***")
    for r, L in lengths.items():
        print(f"    Rep {r}: {L} frames ({L * DT_NS:.1f} ns)")
    print("    Averaged plots are truncated to the SHORTEST replica.")
    print("    Drop unfinished replicas from REPS, or wait for them to finish.\n")

n = min(lengths.values())
S = np.stack([rmsd[r][:n] for r in REPS])                   # (rep, frame, part)
t, m, sd = np.arange(n) * DT_NS, S.mean(0), S.std(0)

print("\nAcross replicas (mean +- SD of per-replica means | overall min / max):")
for j, part in enumerate(PARTS):
    a  = S[:, 1:, j]
    pm = a.mean(1)
    print(f"  {part:<9} mean {pm.mean():5.2f} +- {pm.std():.2f} A | "
          f"min {a.min():5.2f} | max {a.max():5.2f} A")

L = S[:, 1:, 2]
print("\n--- LIGAND RMSD ---")
for rep in REPS:
    c = rmsd[rep][1:, 2]
    print(f"  Rep {rep}: mean {c.mean():5.2f} +- {c.std():.2f} | "
          f"min {c.min():5.2f} | max {c.max():5.2f} A")
print(f"  ALL:   mean {L.mean(1).mean():5.2f} +- {L.mean(1).std():.2f} | "
      f"min {L.min():5.2f} | max {L.max():5.2f} A")


def finish(fig, axes, name, ylabel="RMSD (A)"):
    for ax in np.atleast_1d(axes).ravel():
        ax.set_xlabel("Time (ns)", fontsize=13)
        ax.set_ylabel(ylabel, fontsize=13)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{GRAPH}/{TAG}_{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {GRAPH}/{TAG}_{name}.png")


# 1. mean +- SD, backbone / protein / ligand on one panel
fig, ax = plt.subplots(figsize=(14, 5))
for j, (part, c) in enumerate(PARTS.items()):
    pm = S[:, 1:, j].mean(1)
    ax.plot(t, m[:, j], color=c, lw=1.1,
            label=f"{part}: {pm.mean():.2f} +- {pm.std():.2f} A "
                  f"(max {S[:, 1:, j].max():.2f})")
    ax.fill_between(t, m[:, j] - sd[:, j], m[:, j] + sd[:, j], color=c, alpha=0.15)
ax.set_title(f"{NAME} — mean RMSD across {len(REPS)} replicas", fontsize=14)
finish(fig, ax, "rmsd_mean_bb_prot_lig")

# 2. per-replica traces, protein and ligand side by side
fig, axes = plt.subplots(1, 2, figsize=(18, 5))
for ax, (part, j) in zip(axes, [("Protein", 1), ("Ligand", 2)]):
    for rep in REPS:
        y = rmsd[rep][:n, j]
        ax.plot(t, y, color=REP_C[rep], lw=0.8, alpha=0.65,
                label=f"Rep {rep}: {y[1:].mean():.2f} A (max {y[1:].max():.2f})")
    ax.plot(t, m[:, j], color=PROT_COLOR, lw=1.8, label="Mean")
    ax.fill_between(t, m[:, j] - sd[:, j], m[:, j] + sd[:, j],
                    color=PROT_COLOR, alpha=0.18)
    ax.set_title(f"{part} RMSD — per replica", fontsize=14)
finish(fig, axes, "rmsd_per_replica_prot_lig")

pd.DataFrame(np.column_stack([t, m, sd]),
             columns=["time_ns", "bb", "prot", "lig", "bb_sd", "prot_sd", "lig_sd"]
             ).to_csv(f"{GRAPH}/{TAG}_rmsd_mean.csv", index=False)

# ═══════════════════════════ 2. MM-GBSA (OBC2, 0.15 M) ═════════════════════
if RUN_MMGBSA:
    import mdtraj as md
    import openmm
    from openmm import app, unit

    complex_ = f"{BASE}/complex_stripped.prmtop"
    receptor = f"{BASE}/receptor.prmtop"
    ligand   = f"{BASE}/ligand.prmtop"
    for p in (complex_, receptor, ligand):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p} missing — MM-GBSA needs stripped/receptor/ligand topologies")

    _names = [openmm.Platform.getPlatform(i).getName()
              for i in range(openmm.Platform.getNumPlatforms())]
    if "CUDA" in _names:
        PLATFORM, PROPS, N_SAMPLE = "CUDA", {"Precision": "double"}, N_SAMPLE_GPU
    else:
        PLATFORM, PROPS, N_SAMPLE = "CPU", {"Threads": str(os.cpu_count())}, N_SAMPLE_CPU

    def make_context(prmtop_path):
        system = app.AmberPrmtopFile(prmtop_path).createSystem(
            nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2,
            implicitSolventSaltConc=0.15 * unit.moles / unit.liter)
        integ = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        return openmm.Context(system, integ,
                              openmm.Platform.getPlatformByName(PLATFORM), PROPS), integ

    def gb_energy(ctx, pos):
        ctx.setPositions(pos)
        return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilocalories_per_mole)

    n_complex  = app.AmberPrmtopFile(complex_).topology.getNumAtoms()
    n_receptor = app.AmberPrmtopFile(receptor).topology.getNumAtoms()
    ctx = None

    results, times = [], []
    for rep in REPS:
        legs, imaged = inputs[rep]
        f = f"{CACHE}/mmgbsa_rep{rep}{sfx(imaged)}.npz"
        if os.path.exists(f):
            d = np.load(f)
            results.append(d["dG"]); times.append(d["t"])
            print(f"Rep {rep}: cache, mean dG = {d['dG'].mean():.2f} "
                  f"+- {d['dG'].std():.2f} kcal/mol")
            continue
        if ctx is None:                      # build contexts only if something needs computing
            print(f"MM-GBSA platform: {PLATFORM}, {N_SAMPLE} frames/replica")
            ctx = [make_context(p) for p in (complex_, receptor, ligand)]
        (ctx_c, _), (ctx_r, _), (ctx_l, _) = ctx

        n_tot  = mda.Universe(TOP, legs).trajectory.n_frames
        stride = max(1, n_tot // N_SAMPLE)
        traj   = md.load(legs, top=TOP, stride=stride)
        if not imaged:
            traj.image_molecules(inplace=True)
        traj = traj.atom_slice(range(n_complex))

        dG = np.empty(traj.n_frames)
        for j in range(traj.n_frames):
            pos = traj.xyz[j]
            dG[j] = (gb_energy(ctx_c, pos)
                     - gb_energy(ctx_r, pos[:n_receptor])
                     - gb_energy(ctx_l, pos[n_receptor:n_complex]))
        results.append(dG)
        times.append(np.arange(traj.n_frames) * stride * DT_NS)
        np.savez(f, dG=dG, t=times[-1])
        print(f"Rep {rep}: computed, {len(dG)} frames (stride {stride}), "
              f"mean dG = {dG.mean():.2f} +- {dG.std():.2f} kcal/mol")
    del ctx

    mg      = min(map(len, results))
    G       = np.vstack([x[:mg] for x in results])
    g_m, g_sd, g_t = G.mean(0), G.std(0), times[0][:mg]

    fig, ax = plt.subplots(figsize=(14, 6))
    for rep, dG, tt in zip(REPS, results, times):
        ax.plot(tt, dG, color=REP_C[rep], lw=0.8,
                label=f"Rep {rep}: {dG.mean():.1f} kcal/mol")
    ax.set_title(f"{NAME} — MM-GBSA per replica", fontsize=14)
    finish(fig, ax, "mmgbsa_per_replica", ylabel="dG MM-GBSA (kcal/mol)")

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(g_t, g_m, color=PROT_COLOR, lw=1.5,
            label=f"{NAME}: {G.mean():.1f} +- {G.mean(1).std():.1f} kcal/mol")
    ax.fill_between(g_t, g_m - g_sd, g_m + g_sd, color=PROT_COLOR, alpha=0.2)
    ax.set_title(f"{NAME} — MM-GBSA, mean +- SD across {len(results)} replicas",
                 fontsize=14)
    finish(fig, ax, "mmgbsa_mean", ylabel="dG MM-GBSA (kcal/mol)")

    pd.DataFrame({"time_ns": g_t, "dG_mean": g_m, "dG_sd": g_sd,
                  **{f"rep{r}": G[i] for i, r in enumerate(REPS)}}
                 ).to_csv(f"{GRAPH}/{TAG}_mmgbsa.csv", index=False)

    print(f"\nOverall MM-GBSA: {G.mean():.2f} +- {G.mean(1).std():.2f} kcal/mol "
          f"(SD across replica means)")

print(f"\nFigures and CSVs in: {GRAPH}")
