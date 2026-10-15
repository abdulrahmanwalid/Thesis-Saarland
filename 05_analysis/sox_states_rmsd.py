# Per-residue RMSD comparison of KSHV SOX across three crystal states
#   Apo / DNA-free : 3FHD
#   DNA-bound      : 3POV
#   RNA-bound      : 5HSW
# Step 1: align each pair (chain A) in PyMOL and compute per-residue RMSD of matched atoms
# Step 2: print all per-residue RMSDs
# Step 3: find 3POV residues near DNA and report their RMSDs
# Step 4: scatter plots of residue RMSDs between comparisons (Pearson r, p, y = x line)
# Requires PyMOL as a Python module (pymol-open-source).

# %% Warnings
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
warnings.filterwarnings("ignore")

# %% Download the PDB files from RCSB (skipped if already downloaded)
import os, urllib.request

PDB_DIR = "pdb_files"
os.makedirs(PDB_DIR, exist_ok=True)
for pdb_id in ["3fhd", "3pov", "5hsw"]:
    path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    if not os.path.exists(path):
        print(f"Downloading {pdb_id.upper()} from RCSB...")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", path)
    print(f"{pdb_id.upper()}: {path}")

# %% Step 1: overall and per-residue RMSD
from pymol import cmd
import math

def step_1_average_rmsd():
    """
    Loads structures, calculates overall average RMSD per pair.
    Returns: all_pair_rmsd_dicts — { (name1, name2): { res_label: rmsd } }
             Used by Steps 2, 3, and 4.
    """
    cmd.reinitialize()

    for pdb_id in ["3fhd", "3pov", "5hsw"]:
        print(f"Loading {pdb_id} into memory...")
        cmd.load(os.path.join(PDB_DIR, f"{pdb_id}.pdb"), pdb_id)


    comparison_pairs = [("3fhd", "3pov"), ("3pov", "5hsw"), ("3fhd", "5hsw")]
    all_pair_rmsd_dicts = {}

    print(f"\n{'='*55}")
    print("STEP 1: OVERALL AVERAGE RMSD AND MATCHED LENGTHS")
    print(f"{'='*55}")

    for name1, name2 in comparison_pairs:
        mapping_object = f"map_{name1}_{name2}"
        cmd.align(f"{name1} and chain A", f"{name2} and chain A",
                  cycles=0, object=mapping_object)
        raw_aln = cmd.get_raw_alignment(mapping_object)

        residue_storage = {}
        for pair in raw_aln:
            m1 = cmd.get_model(f"{pair[0][0]} and index {pair[0][1]}")
            m2 = cmd.get_model(f"{pair[1][0]} and index {pair[1][1]}")
            a1, a2 = m1.atom[0], m2.atom[0]
            res_label = f"{a1.resn}{a1.resi}"
            if res_label not in residue_storage:
                residue_storage[res_label] = []
            residue_storage[res_label].append((a1, a2))

        rmsd_dict = {}
        total_rmsd_sum = 0.0
        for res_label, atom_list in residue_storage.items():
            sos = sum(
                (a1.coord[0]-a2.coord[0])**2 +
                (a1.coord[1]-a2.coord[1])**2 +
                (a1.coord[2]-a2.coord[2])**2
                for a1, a2 in atom_list
            )
            rmsd = math.sqrt(sos / len(atom_list))
            rmsd_dict[res_label] = rmsd
            total_rmsd_sum += rmsd

        matched_count = len(rmsd_dict)
        avg_rmsd = total_rmsd_sum / matched_count if matched_count > 0 else 0.0
        all_pair_rmsd_dicts[(name1, name2)] = rmsd_dict

        print(f"{name1.upper()} vs {name2.upper()}:")
        print(f"  - Total Matched Residues (Length): {matched_count}")
        print(f"  - Overall Average RMSD: {avg_rmsd:.3f} Å\n")

    return all_pair_rmsd_dicts

all_pair_rmsd_dicts = step_1_average_rmsd()


# %% Step 2: print per-residue RMSDs
def step_2_individual_rmsd(all_pair_rmsd_dicts):
    """
    Prints sorted per-residue RMSDs for all three pairs.
    Uses the pre-computed dicts from Step 1 — no PyMOL calls needed.
    """
    pairs_to_compare = [("3fhd", "3pov"), ("3fhd", "5hsw"), ("3pov", "5hsw")]

    for name1, name2 in pairs_to_compare:
        rmsd_dict = all_pair_rmsd_dicts[(name1, name2)]

        print(f"\n{'='*60}")
        print(f"STEP 2: ALL INDIVIDUAL RMSD FOR {name1.upper()} vs {name2.upper()}")
        print(f"{'='*60}")
        print(f"{'Residue':<12} {'RMSD (Å)':>10}")
        print("-" * 40)

        # Sort numerically by residue number extracted from label
        def sort_key(label):
            digits = ''.join(filter(str.isdigit, label))
            return int(digits) if digits else 0

        for res_label in sorted(rmsd_dict, key=sort_key):
            print(f"{res_label:<12} {rmsd_dict[res_label]:>10.2f}")

        print("-" * 40)
        print(f"Total matched residues (length) for {name1} vs {name2}: {len(rmsd_dict)}")

step_2_individual_rmsd(all_pair_rmsd_dicts)


# %% Step 3: DNA-contact residues
def step_3_dna_target_rmsd(all_pair_rmsd_dicts, distance_cutoff=5.0):
    literature_sites = {"ASP221", "GLU244", "LYS246", "GLU184"}

    print(f"\n{'='*60}")
    print(f"EXTRACTING DNA INTERACTING SITES (Cutoff: {distance_cutoff} Å)")
    print(f"{'='*60}")

    # DNA still in memory here — do the selection FIRST
    cmd.select("bindC", f"byres 3pov and chain A within {distance_cutoff} of 3pov and chain C")
    cmd.select("bindD", f"byres 3pov and chain A within {distance_cutoff} of 3pov and chain D")
    cmd.select("dna_contacts", "bindC or bindD and name CA")

    dynamic_dna_residues = []
    cmd.iterate("dna_contacts", "dynamic_dna_residues.append(resn+resi)",
                space={"dynamic_dna_residues": dynamic_dna_residues})
    dynamic_dna_residues = list(set(dynamic_dna_residues))
    target_residues = [r for r in dynamic_dna_residues if r not in literature_sites]

    print(f"Found {len(dynamic_dna_residues)} total DNA-interacting residues.")
    print(f"Target residues after excluding literature sites: {len(target_residues)}")
    print(target_residues)

    # NOW clean for alignment
    cmd.remove("not polymer.protein")

    pairs_to_compare = [("3fhd", "3pov"), ("3fhd", "5hsw"), ("3pov", "5hsw")]

    for name1, name2 in pairs_to_compare:
        rmsd_dict = all_pair_rmsd_dicts[(name1, name2)]
        matched = {r: rmsd_dict[r] for r in target_residues if r in rmsd_dict}

        print(f"\n{'='*60}")
        print(f"TARGET RESIDUES RMSD: {name1.upper()} vs {name2.upper()}")
        print(f"{'='*60}")
        print(f"{'Residue':<12} {'RMSD (Å)':>10}")
        print("-" * 40)

        def sort_key(label):
            digits = ''.join(filter(str.isdigit, label))
            return int(digits) if digits else 0

        for res_label in sorted(matched, key=sort_key):
            print(f"{res_label:<12} {matched[res_label]:>10.2f}")

        print("-" * 40)
        avg = sum(matched.values()) / len(matched) if matched else 0.0
        print(f"Total matched target residues: {len(matched)}")
        print(f"Average RMSD for these targets: {avg:.2f} Å")

    return target_residues

target_residues = step_3_dna_target_rmsd(all_pair_rmsd_dicts, distance_cutoff=5.0)


# %% Step 4: scatter plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

def step_4_scatter_plot(all_pair_rmsd_dicts, target_residues, out_png="sox_states_rmsd_scatter.png"):
    """
    Scatter plot using pre-computed data from Steps 1 & 3. No PyMOL calls.
    Lit sites and DNA-contact residues are labelled; the dashed line is y = x.
    """
    literature_sites = {"ASP221", "GLU244", "LYS246", "GLU184"}

    rdf = all_pair_rmsd_dicts[("3fhd", "3pov")]   # DNA-bound vs free
    rrf = all_pair_rmsd_dicts[("3fhd", "5hsw")]   # RNA-bound vs free
    rdr = all_pair_rmsd_dicts[("3pov", "5hsw")]   # DNA-bound vs RNA-bound

    pymol_nums = set(target_residues)
    lit_nums   = set(literature_sites)
    common     = list(set(rdf) & set(rrf) & set(rdr))

    lit_res   = [r for r in common if r in lit_nums]
    pymol_res = [r for r in common if r in pymol_nums]
    norm_res  = [r for r in common if r not in lit_nums and r not in pymol_nums]

    def get_xy(res_list, dx, dy):
        return [dx[r] for r in res_list], [dy[r] for r in res_list]

    comparisons = [
        (rdf, rdr, "DNA Free", "DNA·RNA"),
        (rrf, rdr, "RNA Free", "DNA·RNA"),
        (rdf, rrf, "DNA Free", "RNA Free"),
    ]

    # same axis range on all panels so the y = x line is comparable
    lim = max(max(d[r] for r in common) for d in (rdf, rrf, rdr)) * 1.08

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    for ax, (dx, dy, xlab, ylab) in zip(axes, comparisons):
        if norm_res:
            x, y = get_xy(norm_res, dx, dy)
            ax.scatter(x, y, c="lightgray", s=30, marker="o", zorder=1)
        if pymol_res:
            x, y = get_xy(pymol_res, dx, dy)
            ax.scatter(x, y, c="darkorange", s=60, marker="s", zorder=3)
            for r, xi, yi in zip(pymol_res, x, y):
                ax.annotate(r, (xi, yi), xytext=(5, 4), textcoords="offset points", fontsize=7,
                            bbox=dict(boxstyle="round,pad=0.15", fc="lightyellow", ec="gray", lw=0.5),
                            zorder=3)
        if lit_res:
            x, y = get_xy(lit_res, dx, dy)
            ax.scatter(x, y, c="firebrick", s=220, marker="*", edgecolors="k", lw=0.5, zorder=4)
            for r, xi, yi in zip(lit_res, x, y):
                ax.annotate(r, (xi, yi), xytext=(6, -8), textcoords="offset points",
                            fontsize=9, fontweight="bold",
                            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="k"), zorder=5)

        # perfect correlation line (y = x)
        ax.plot([0, lim * 0.96], [0, lim * 0.96], "k--", lw=1.5, alpha=0.8, zorder=2)

        all_x, all_y = get_xy(common, dx, dy)
        r_val, p_val = pearsonr(all_x, all_y)
        ax.text(0.02, 0.98, f"r = {r_val:.3f}\np = {p_val:.2e}",
                transform=ax.transAxes, va="top",
                bbox=dict(facecolor="white", edgecolor="k", alpha=0.9))

        ax.set_xlim(-0.3, lim)
        ax.set_ylim(-0.3, lim)
        ax.set_xlabel(f"RMSD {xlab} (Å)")
        ax.set_ylabel(f"RMSD {ylab} (Å)")
        ax.grid(True, ls="--", alpha=0.5)

    handles = [
        plt.Line2D([], [], marker="o", ls="", color="lightgray", markersize=8, label="Normal"),
        plt.Line2D([], [], marker="*", ls="", color="firebrick", markersize=14, label="Lit Site"),
        plt.Line2D([], [], marker="s", ls="", color="darkorange", markersize=8, label="PyMOL DNA"),
        plt.Line2D([], [], ls="--", color="k", label="Perfect correlation (y=x)"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=11, frameon=True)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()

step_4_scatter_plot(all_pair_rmsd_dicts, target_residues)