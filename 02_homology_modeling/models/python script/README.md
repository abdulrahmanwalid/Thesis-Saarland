# Homology modeling

Homology models of three herpesvirus alkaline nucleases built in MODELLER on the KSHV SOX structure (3POV) as template. The template is the fully loop-modelled 3POV used in the MD simulations, so modelled loops carry over into the homologs. The two catalytic Mg2+ ions and their coordinating waters are transferred into each model.

## Files

| File | Purpose |
|---|---|
| `clean_template.pml` | PyMOL script: strips DNA, hydrogens and salt from the MD 3POV, keeps protein + 2 Mg2+ + waters within 2.8 Å of Mg, saves `3pov_tmpl.pdb` |
| `modeller_homolog.py` | Fetches the target sequence from UniProt, aligns it to the template (`align2d`) and builds the models (`AutoModel`) |
| `models/` | Final selected model for each homolog |

## Targets

| Homolog | Virus | UniProt | Length (aa) | Catalytic residues |
|---|---|---|---|---|
| UL12 | HSV-1 | G8HBC5 | 626 | E280, D340, E364, K366 |
| UL98 | HCMV | F5HF49 | 584 | D254, E278, K280 |
| BGLF5 | EBV | P03217 | 470 | E166, D203, E225, K227 |

Template catalytic residues (KSHV SOX, 3POV): E184, D221, E244, K246.

## How to run

1. Open `clean_template.pml`, set `<YOUR_3POV>` to your loop-modelled 3POV file, and run it in PyMOL (`@clean_template.pml`)
2. In `modeller_homolog.py`, fill in `UNIPROT`, `NAME` and `TMPL_DIR` from the table above
3. Run it, then open `<NAME>.pap` and check the catalytic residues align to the SOX residues
4. Repeat steps 2–3 for each homolog

## Model selection

Models are ranked by DOPE score, but the final model was not chosen on DOPE alone. Each candidate was inspected visually, and the selected model is the one with the best overall shape: intact fold, sensible loops, and a well-formed active site with the Mg2+ and catalytic residues in place.
