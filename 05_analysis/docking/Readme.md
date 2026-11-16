# Docking analysis — RMSD and MM-GBSA

Trajectory analysis of AD-51 docked into KSHV SOX (PDB 3POV) and, separately,
into three herpesvirus alkaline nuclease homologs. The SOX analysis compares
four ligand states at two sites; the homolog analysis compares four proteins at
one ligand state. The two are kept in separate folders because they answer
different questions and use different figure colours.

## Pipeline

1. `make_mmgbsa_topologies.in` — cpptraj, once per system, builds the three
   solute topologies MM-GBSA needs from the solvated complex topology.
2. `rmsd_mmgbsa.py` — once per system. Edit the CONFIG block at the top, run it,
   and it writes per-replica results to that system's `cache/`.
3. `sox/compiled_rmsd.py` and `sox/compiled_mmgbsa.py` — read those caches and
   produce the figures. They compute nothing themselves.

Trajectories and topologies are not included here due to size.

## Method

Production MD on an HPC cluster, 5 replicas per system, 200 ns each, 10 ps per
frame. The first 1.12 ns of each trajectory is equilibration and is excluded
from every reported average; the caches hold raw per-frame values, so the cutoff
is applied at the figure stage.

**RMSD** — protein backbone, whole protein, and ligand, fitted on the protein
backbone with frame 0 as reference. Frame 0 is excluded from reported means,
since its RMSD is zero by definition.

**MM-GBSA** — OBC2 implicit solvent at 0.15 M salt, single-trajectory
approximation: dG = E(complex) - E(receptor) - E(ligand), all three evaluated on
the same frames. Mg2+ is kept with the receptor, so the result is the binding
free energy of AD-51 to the metal-bound enzyme rather than to the apo protein.

Because the ligand leaves its docked site in some systems, the compiled MM-GBSA
script also reports a value restricted to frames where the ligand stays within
3 A of its docked position, alongside the fraction of frames that satisfies it.
These figures measure association with the complex rather than occupancy of a
defined site, and should be read accordingly.

Replicas are averaged after truncation to the shortest; `rmsd_mmgbsa.py` warns
if replica lengths differ by more than 2%. Trajectories needing periodic-boundary
correction were passed through `cpptraj autoimage` first — the scripts pick up
`trajectory_imaged.dcd` automatically and check its frame count against the raw
trajectory.

## SOX systems

| System | BASE | Replica pattern | Reps |
|---|---|---|---|
| Thr189 site, neutral lactam | `3POV/Thr189_dock` | `Thr189_Dock_{r}_3POV` | 5 |
| Mg2+ site, neutral lactam | `3POV/Mg_dock` | `Mg_Dock_{r}_3POV` | 5 |
| Mg2+ site, lactim tautomer | `3POV/Mg_dock/Mg_Docked_Poses/Neutral_Charge_Pose` | `Mg_Dock_{r}_Neutral_Pose_3POV` | 5 |
| Mg2+ site, deprotonated | `3POV/Mg_dock/Mg_Docked_Poses/Negative_Charge_Pose` | `Mg_Dock_{r}_Negative_Pose` | 5 |

All four have MM-GBSA. The lactam and lactim are the two neutral tautomers; the
deprotonated form carries its negative charge delocalised over the amide N-C=O
unit rather than on a single atom.

Folder names on disk use the earlier neutral/negative wording, which predates the
tautomer assignment — the table above is authoritative for which is which.

## Homolog systems

| System | BASE | Replica pattern | Reps |
|---|---|---|---|
| EBV BGLF5 | `EBV/EBV_Docking/AD51 Negative Nitrogen` | `Replica{r}` | 5 |
| HSV-1 UL12 | `UL12/UL12 Docking/AD51 Negative Nitrogen` | `Replica{r}` | 5 |
| HCMV UL98 | `UL98/UL98 Docking/UL98 Negative Nitrogen` | `Replica{r}` | 5 |

All three use the deprotonated ligand at the Mg2+ site. RMSD only; MM-GBSA was
not run for the homologs.

## Figure colours

The two figure families use the same four colours for different things, so read
each within its own section.

SOX figures colour by **ligand state**:

| State | Colour |
|---|---|
| Thr189 site, neutral lactam | `tab:blue` |
| Mg2+ site, neutral lactam | `tab:orange` |
| Mg2+ site, lactim tautomer | `tab:green` |
| Mg2+ site, deprotonated | `tab:red` |

Homolog figures colour by **protein**: SOX `tab:blue`, EBV `tab:orange`,
UL12 `tab:green`, UL98 `tab:red`. Backbone traces are `steelblue` and ligand
traces `crimson` throughout. Per-replica traces use a separate neutral palette
generated from the replica count, so it cannot fall out of step when a replica
is added.

## Software

MDAnalysis, MDTraj, OpenMM, AmberTools (cpptraj), NumPy, pandas, Matplotlib.
