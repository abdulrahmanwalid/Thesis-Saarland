# DCCM: Baseline vs PTR253

Dynamic cross-correlation analysis (DCCM) of KSHV SOX, comparing the unphosphorylated system (baseline) with the phosphotyrosine 253 system (PTR). Both are built from 3POV with DNA and the two catalytic Mg2+ ions, 3 replicas each.

## Notebook

`dccm_ptr.ipynb` contains the code and the plots:

1. **Heatmaps**: replica-averaged DCCM for baseline and PTR, and the difference (PTR − Baseline), with residue 253 marked
2. **Chord diagrams**: strong correlations in each system, and the changes between them drawn as chords around the sequence

## Method

- **Alignment**: trajectories are aligned on the same stable core as the RMSD/RMSF analysis (30% lowest-RMSF residues of baseline replica 1)
- **DCCM**: Cα fluctuations around the mean structure; C_ij = ⟨Δr_i·Δr_j⟩ / √(⟨Δr_i²⟩⟨Δr_j²⟩), averaged over the 3 replicas
- **Chords (per system)**: pairs with r ≥ 0.5 (red) or r ≤ −0.5 (blue), at least 4 residues apart
- **Chords (difference)**: only robust changes, where |Δ| is larger than the summed replica SDs and at least 0.35; red = strengthened, blue = weakened
- **Numbering**: residue numbers follow the topology, which is +2 relative to the crystal structure

## Requirements

MDAnalysis, NumPy, Matplotlib. The notebook was run on Google Colab with the data on Google Drive; trajectories (`.dcd`) and topologies (`.prmtop`) are not included in the repository, so update `BASE` to point to your own files.
