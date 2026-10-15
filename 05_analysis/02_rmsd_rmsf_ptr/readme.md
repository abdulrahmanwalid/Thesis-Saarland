# RMSD and RMSF: Baseline vs PTR253

Compares the dynamics of unphosphorylated KSHV SOX (baseline) with the phosphotyrosine 253 system (PTR), both built from 3POV with DNA and the two catalytic Mg2+ ions. Each system was simulated in 3 independent replicas.

## Notebook

`rmsd_rmsf_ptr.ipynb` contains the code and the plots:

1. **RMSD per region**: mean ± SD over the 3 replicas for
   - the full protein (Cα)
   - the catalytic site: E186, D223, E246, K248 and the Mg2+ ions
   - the PTR253 region (residues 284–289)
2. **Per-residue RMSF**: Cα RMSF for both systems, with the catalytic residues and the PTR region marked, plus a ΔRMSF panel (PTR − Baseline)

## Method

- **Stable-core alignment**: the 30% of residues with the lowest RMSF in baseline replica 1 define a stable core. Both systems are aligned on this same core, so flexible loops do not distort the comparison.
- **RMSD**: each replica is aligned on the core and compared to its first frame.
- **RMSF**: two-pass alignment, first on the core, then on the average structure.
- **Numbering**: residue numbers follow the topology, which is +2 relative to the crystal structure (e.g. E186 here = E184 in 3POV).

## Requirements

MDAnalysis, NumPy, Matplotlib. Trajectories (`.dcd`) and topologies (`.prmtop`) are not included in the repository; update the paths at the top of each cell to point to your own files.
