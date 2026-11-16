# Docking analysis — RMSD and MM-GBSA

Trajectory analysis of AD51 docked into KSHV SOX and three herpesvirus alkaline
nuclease homologs. All systems were analysed with the same script,
[`rmsd_mmgbsa.py`](rmsd_mmgbsa.py); only the CONFIG block at the top changes
between them. Each system folder holds its own config block and figures.

## Protocols

Production MD was run on an HPC cluster in Saarland Campus, 5 replicas per system, 10 ps per frame.
For each replica the script computes:

- **RMSD** of the protein backbone, the whole protein, and the ligand, fitted on
  the protein backbone with frame 0 as reference. Frame 0 is excluded from all
  reported means, since its RMSD is zero by definition.
- **MM-GBSA** binding free energy in OBC2 implicit solvent at 0.15 M salt, using
  the single-trajectory approximation: the complex, receptor and ligand energies
  are evaluated on the same frames, and dG = E(complex) - E(receptor) - E(ligand).

Replicas are averaged after truncation to the shortest one; the script warns if
replica lengths differ by more than 2%. Per-replica results are cached as `.npy`
and `.npz` under each system's `cache/`, so reruns only recompute what is missing.
Trajectories that needed periodic-boundary correction were passed through
`cpptraj autoimage` first; the script picks up `trajectory_imaged.dcd`
automatically and checks its frame count against the raw trajectory.

Trajectories and topologies are not included in this repository due to size.

## Systems

| Folder | System | Replica pattern | Reps | MM-GBSA |
|---|---|---|---|---|
| `sox_mg_pose1/` | SOX, Mg site, pose 1 | `Mg_Dock_{r}_3POV` | 5 | yes |
| `sox_mg_pose2_neutral/` | SOX, Mg site, pose 2, neutral tautomer | `Mg_Dock_{r}_Neutral_Pose_3POV` | 5 | yes |
| `sox_mg_pose3_negative/` | SOX, Mg site, pose 3, anionic | `Mg_Dock_{r}_Negative_Pose` | 5 | yes |
| `sox_thr189_pose1/` | SOX, Thr189 site, pose 1 | `Thr189_Dock_{r}_3POV` | 5 | yes |
| `sox_thr189_pose2_negative/` | SOX, Thr189 site, pose 2, anionic | `Replica{r}` | 5 | TBD |
| `ebv/` | EBV BGLF5, Mg site, anionic | `Replica{r}` | 5 | no |
| `ul12/` | HSV-1 UL12, Mg site, anionic | `Replica{r}` | 5 | no |
| `ul98/` | HCMV UL98, Mg site, anionic | `Replica{r}` | 5 | no |

The anionic form of AD51 carries its negative charge delocalised over the amide
N-C=O unit rather than on a single atom; the neutral tautomer carries a hydrogen
on the carbonyl oxygen.

## Figure colours

One palette across every figure in the chapter, so a protein keeps its colour
wherever it appears:

| Element | Colour |
|---|---|
| SOX | `tab:blue` |
| EBV BGLF5 | `tab:orange` |
| HSV-1 UL12 | `tab:green` |
| HCMV UL98 | `tab:red` |
| Backbone | `steelblue` |
| Ligand | `crimson` |

Per-replica traces use a separate neutral palette, generated to match the number
of replicas so it cannot fall out of step when a replica is added.

## Running

Edit the CONFIG block and run:

```bash
python rmsd_mmgbsa.py
```

Figures and CSVs are written to `<BASE>/graphs/`, caches to `<BASE>/cache/`.
To force a recompute, delete the relevant cache files.

MM-GBSA needs three extra topologies alongside `complex.prmtop`:
`complex_stripped.prmtop`, `receptor.prmtop` and `ligand.prmtop`, with receptor
atoms ordered before ligand atoms in the stripped topology. Set
`RUN_MMGBSA = False` for RMSD only.

## Software

MDAnalysis, MDTraj, OpenMM, AmberTools (cpptraj), NumPy, pandas, Matplotlib.
