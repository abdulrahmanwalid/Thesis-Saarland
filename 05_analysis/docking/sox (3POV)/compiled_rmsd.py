def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ═══════════════════════════════════════════════════════════════════════════
# Compiled RMSD figure for the four AD-51 complexes.
# Reads the per-replica caches written by rmsd_mmgbsa.py — computes nothing.
# Cache column order: [backbone, protein, ligand]
# ═══════════════════════════════════════════════════════════════════════════

BASE = "/home/abel00002/3POV"
OUT  = f"{BASE}/unified_graphs/rmsd_ad51_complexes.png"
DT, EQ = 0.01, 1.12                 # ns per frame, end of equilibration (ns)
TAIL   = 50.0                       # ns window for the "final 50 ns" summary

SYSTEMS = [
    ("Thr189 site, neutral lactam",
     f"{BASE}/Thr189_dock/cache/rmsd_rep{{k}}.npy", "tab:blue"),
    ("Mg$^{2+}$ site, neutral lactam",
     f"{BASE}/Mg_dock/cache/rmsd_rep{{k}}.npy", "tab:orange"),
    ("Mg$^{2+}$ site, lactim tautomer",
     f"{BASE}/Mg_dock/Mg_Docked_Poses/cache/rmsd_neutral_rep{{k}}.npy", "tab:green"),
    ("Mg$^{2+}$ site, deprotonated",
     f"{BASE}/Mg_dock/Mg_Docked_Poses/cache/rmsd_negative_rep{{k}}.npy", "tab:red"),
]

PANELS = [("Protein backbone", 0), ("Ligand", 2)]


def pick(p, ext):
    im = p.replace(ext, f"_imaged{ext}")
    return im if os.path.isfile(im) else p


print(f"{'system':<34s}{'backbone':>20s}{'ligand':>20s}{'ligand, final 50 ns':>24s}")
print("-" * 98)

results = {}
for label, rpat, color in SYSTEMS:
    traces, bb_means, lig_means, tail_means = [], [], [], []

    for k in range(1, 6):
        r = np.load(pick(rpat.format(k=k), ".npy"))
        t = np.arange(len(r)) * DT
        prod = t >= EQ
        tail = t >= (t.max() - TAIL)

        bb_means.append(r[prod, 0].mean())
        lig_means.append(r[prod, 2].mean())
        tail_means.append(r[tail, 2].mean())
        traces.append((t, r))

    print(f"{label[:34]:<34s}"
          f"{np.mean(bb_means):12.2f} +/- {np.std(bb_means):4.2f}"
          f"{np.mean(lig_means):12.2f} +/- {np.std(lig_means):4.2f}"
          f"{np.mean(tail_means):16.2f} +/- {np.std(tail_means):4.2f}")
    print(f"{'':34s}  replica means, final 50 ns: "
          f"{', '.join(f'{v:.2f}' for v in sorted(tail_means))}")

    results[label] = (traces, color)

print("\nAngstrom, mean across replicas +/- SD between them; "
      f"equilibration (first {EQ} ns) excluded")

# ──────────────────────────────────────────────────────────── figure ───────
fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))

for ax, (panel, col) in zip(axes, PANELS):
    for label, (traces, color) in results.items():
        n = min(len(r) for _, r in traces)
        t = traces[0][0][:n]
        y = np.stack([r[:n, col] for _, r in traces])
        m, s = y.mean(0), y.std(0)
        ax.plot(t, m, color=color, lw=1.2, label=label)
        ax.fill_between(t, m - s, m + s, color=color, alpha=0.18, lw=0)

    ax.axvline(EQ, color="grey", ls=":", lw=0.9)
    ax.set_title(panel)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r"RMSD ($\AA$)")
    ax.set_xlim(0, t.max())
    ax.grid(alpha=0.25, lw=0.6)

axes[0].legend(fontsize=9, frameon=False, loc="lower right")

fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"\nsaved {OUT}")
