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

BASE = "/home/abel00002/3POV"
OUT  = f"{BASE}/unified_graphs/mmgbsa_ad51_complexes.png"
DT, EQ, BOUND = 0.01, 1.12, 3.0     # ns/frame, equilibration, ligand-RMSD cutoff

SYSTEMS = [
    ("Thr189 site, neutral lactam",
     f"{BASE}/Thr189_dock/cache/mmgbsa_rep{{k}}.npz",
     f"{BASE}/Thr189_dock/cache/rmsd_rep{{k}}.npy", "tab:blue"),
    ("Mg$^{2+}$ site, neutral lactam",
     f"{BASE}/Mg_dock/cache/mmgbsa_rep{{k}}.npz",
     f"{BASE}/Mg_dock/cache/rmsd_rep{{k}}.npy", "tab:orange"),
    ("Mg$^{2+}$ site, lactim tautomer",
     f"{BASE}/Mg_dock/Mg_Docked_Poses/cache/mmgbsa_neutral_rep{{k}}.npz",
     f"{BASE}/Mg_dock/Mg_Docked_Poses/cache/rmsd_neutral_rep{{k}}.npy", "tab:green"),
    ("Mg$^{2+}$ site, deprotonated",
     f"{BASE}/Mg_dock/Mg_Docked_Poses/cache/mmgbsa_negative_rep{{k}}.npz",
     f"{BASE}/Mg_dock/Mg_Docked_Poses/cache/rmsd_negative_rep{{k}}.npy", "tab:red"),
]

def pick(p, ext):
    im = p.replace(ext, f"_imaged{ext}")
    return im if os.path.isfile(im) else p

print(f"{'system':<32s}{'full run':>22s}{'bound frames only':>24s}{'bound':>8s}")
print("-" * 88)

results = {}
for label, gpat, rpat, color in SYSTEMS:
    per_rep_full, per_rep_bound, traces, kept, total = [], [], [], 0, 0

    for k in range(1, 6):
        g  = np.load(pick(gpat.format(k=k), ".npz"))
        dG, t = g["dG"], g["t"]
        rms = np.load(pick(rpat.format(k=k), ".npy"))[:, 2]      # ligand column

        idx   = np.clip(np.round(t / DT).astype(int), 0, len(rms) - 1)
        prod  = t >= EQ
        bound = prod & (rms[idx] < BOUND)

        per_rep_full.append(dG[prod].mean())
        if bound.sum() > 20:
            per_rep_bound.append(dG[bound].mean())
        kept  += bound.sum()
        total += prod.sum()
        traces.append((t, dG))

    f_m, f_s = np.mean(per_rep_full), np.std(per_rep_full)
    if per_rep_bound:
        b = f"{np.mean(per_rep_bound):8.2f} +/- {np.std(per_rep_bound):5.2f}"
    else:
        b = "          n/a       "
    print(f"{label[:32]:<32s}{f_m:14.2f} +/- {f_s:5.2f}{b:>24s}{100*kept/total:7.1f}%")

    results[label] = (per_rep_full, traces, color)

print("\nkcal/mol, mean across replicas +/- SD between them")
print(f"bound frames = ligand RMSD < {BOUND} A; 'bound' column is the share of production frames kept")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 4.8),
                             gridspec_kw={"width_ratios": [2, 1]})

for label, (per_rep, traces, color) in results.items():
    n = min(len(t) for t, _ in traces)
    t = traces[0][0][:n]
    y = np.stack([d[:n] for _, d in traces])
    m, s = y.mean(0), y.std(0)
    a1.plot(t, m, color=color, lw=1.2, label=label)
    a1.fill_between(t, m - s, m + s, color=color, alpha=0.18, lw=0)

labels = list(results)
for i, label in enumerate(labels):
    per_rep, _, color = results[label]
    a2.scatter([i] * len(per_rep), per_rep, color=color, s=28, zorder=3)
    a2.errorbar(i, np.mean(per_rep), yerr=np.std(per_rep), color=color,
                capsize=5, lw=1.4, zorder=2)

a1.axvline(EQ, color="grey", ls=":", lw=0.9)
a1.set_xlabel("Time (ns)")
a1.set_ylabel(r"$\Delta G$ (kcal mol$^{-1}$)")
a1.set_xlim(0, t.max())
a1.legend(fontsize=9, frameon=False)
a1.grid(alpha=0.25, lw=0.6)

a2.set_xticks(range(len(labels)))
a2.set_xticklabels(["Thr189\nlactam", "Mg\nlactam", "Mg\nlactim", "Mg\ndeprot."], fontsize=9)
a2.set_ylabel(r"$\Delta G$ (kcal mol$^{-1}$)")
a2.set_title("Replica means")
a2.grid(alpha=0.25, lw=0.6, axis="y")

fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"\nsaved {OUT}")
