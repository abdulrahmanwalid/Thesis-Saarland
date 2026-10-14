# tleap inputs

Reconstructed tleap inputs for each simulation system (originally typed interactively).
Run with `tleap -f <file>.in`.

Input and output file names (PDB, mol2, frcmod, prmtop) are generic placeholders; replace them with your own.

| File | System |
|---|---|
| `sox_baseline.in` | SOX + DNA + Mg2+ |
| `sox_ptr253.in` | SOX with phosphotyrosine 253 |
| `sox_ad51_mgsite.in` | SOX + AD51 at the Mg2+ site |
| `sox_ad51_mgsite_negative_nitrogen.in` | Mg2+ site pose, AD51 with deprotonated ring nitrogen (N−) |
| `sox_ad51_mgsite_lactim.in` | Mg2+ site pose, AD51 lactim tautomer (OH) |
| `sox_ad51_thr189.in` | SOX + AD51 at the Thr189/Tyr243 site |
| `homologs_dna.in` | UL12, UL98 and BGLF5 homology models + DNA + Mg2+ (one file, swap names) |

Force fields: ff14SB (protein), phosaa14SB (PTR), OL15 (DNA), GAFF2/AM1-BCC (AD51), TIP3P water. The Mg2+ site poses explicitly load the Li/Merz 12-6 Mg2+ parameters (`frcmod.ions234lm_126_tip3p`).

## Ion counts (SPLIT method)

Ions follow the SPLIT method (Machado & Pantano, 2020) for 0.15 M NaCl:

1. Run the file once up to `solvatebox` and note the number of waters (N_w) and the system charge (Q) from `charge`
2. N0 = 0.15 × N_w / 56
3. Na+ = round(N0 − Q/2), Cl− = round(N0 + Q/2)
4. Put these into `addionsrand` and run the full file
