# SOX crystal states comparison

Compares three crystal structures of KSHV SOX to see whether the active site changes shape depending on what is bound:

| PDB | State |
|---|---|
| 3FHD | Apo (no nucleic acid) |
| 3POV | DNA-bound |
| 5HSW | RNA-bound |

## What the script does

1. Downloads the three structures from RCSB
2. Aligns each pair (chain A) in PyMOL and computes the RMSD of every matched residue
3. Finds the 3POV residues within 5 Å of the DNA
4. Plots residue RMSDs from one comparison against another (DNA Free, RNA Free, DNA·RNA), with Pearson r and p and a y = x reference line

Marked in the plot:
- **Lit Site** (red stars): catalytic residues E184, D221, E244, K246
- **PyMOL DNA** (orange squares): DNA-contacting residues in 3POV

Residue numbering follows the crystal structures.

## How to run

```bash
conda create -n pymol_env -c conda-forge python=3.10 pymol-open-source scipy matplotlib -y
conda activate pymol_env
python sox_states_rmsd.py
```

Output: per-residue RMSD tables printed to the console and `sox_states_rmsd_scatter.png`.
