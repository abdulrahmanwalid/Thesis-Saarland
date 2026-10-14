# Computational Analysis of Non-Competitive Inhibitors of Alkaline Nucleases

MSc thesis in Bioinformatics, Saarland University, carried out at the Helmholtz Institute for Pharmaceutical Research Saarland (HIPS).

**Author:** Abdulrahman Walid
**Supervisors:** Prof. Dr. Andrea Volkamer and Prof. Dr. Martin Empting
**Advisor:** Hamza Ibrahim (PhD student)

## About

Herpesviruses encode an alkaline nuclease that is essential for viral replication, which makes it an attractive antiviral target. This thesis uses the KSHV SOX protein (PDB: 3POV) as a structural model to study how the small molecule AD51 can inhibit the enzyme non-competitively, without competing with DNA for the active site.

The work combines:
- **System preparation** of SOX with and without DNA, ligand parameterization
- **Molecular docking** of AD51 at candidate binding sites
- **Molecular dynamics simulations** (AMBER force fields, OpenMM)
- **Binding free energy estimation** (MM-GBSA)
- **Homology modeling** of related herpesvirus nucleases (HSV-1 UL12, HCMV UL98, EBV BGLF5) to compare binding sites across the family

## Repository structure

| Folder | Contents |
|---|---|
| `01_system_preparation/` | Protocols and tleap inputs for building the simulation systems |
| `02_homology_modeling/` | MODELLER workflow and final homolog models |
| `03_docking/` | Docking protocol and selected poses |
| `04_md_simulations/` | Minimization and production MD notebooks, HPC job scripts |
| `05_analysis/` | Trajectory analysis notebooks (RMSD/RMSF, DCCM, contacts, MM-GBSA) |
| `06_visualization/` | VMD and PyMOL scripts |

Trajectories and topologies are not included due to size.

## Tools

AmberTools (tleap, antechamber), OpenMM, MDAnalysis, MODELLER, MOE, PyMOL, VMD
