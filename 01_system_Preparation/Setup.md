# System preparation

Steps taken to turn the 3POV crystal structure into simulation-ready systems. These were done interactively (ChimeraX, PyMOL, MOE, command line), so they are described here rather than scripted. The final tleap inputs are in `tleap_inputs/`.

## 1. Loop modeling

- Removed DNA (chains C/D) and the formate (FMT) crystallization artifact from 3POV in ChimeraX, keeping the protein, the two Mg2+ ions and crystal waters
- Built the missing loops with Modeller (ChimeraX Model Loops), using the RCSB construct sequence
- Chose the final loop model by visual inspection and clash checks, not zDOPE alone
- Saved with Amber numbering (continuous residue numbers for tleap)
- Superimposed the finished model back onto 3POV in PyMOL and copied the DNA back in

## 2. Protonation

- Stripped crystal waters, then added hydrogens with Reduce
- Cleaned the structure with pdb4amber for AMBER residue and atom names
- For the phosphorylated system, Tyr253 was renamed to PTR (parameters from phosaa14SB)
- For docked complexes, receptor hydrogens were removed and pdb4amber was run with `--nohyd`, letting tleap add hydrogens to avoid naming mismatches

## 3. Ligand parameters (AD51)

- Exported the docked AD51 pose from MOE as MOL2 (keeps bond orders and charge state), residue name `LIG`
- Assigned GAFF2 atom types and AM1-BCC charges with antechamber, setting the net charge explicitly:
  - lactim (neutral): `-nc 0`
  - negative nitrogen (N−): `-nc -1`

  ```bash
  antechamber -i lig.mol2 -fi mol2 -o lig_gaff.mol2 -fo mol2 -c bcc -nc <charge> -at gaff2 -rn LIG
  parmchk2 -i lig_gaff.mol2 -f mol2 -o lig.frcmod -s gaff2
  ```
- Checked the total charge of the output and the missing parameters in the frcmod

## 4. Topology

See `tleap_inputs/` for the tleap files and the SPLIT ion calculation.
