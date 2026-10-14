# Clean the loop-modelled 3POV (from the MD setup) into a MODELLER template:
# protein + 2 Mg2+ + Mg-coordinating waters, no DNA.
# Run in PyMOL:  @clean_template.pml

# fill in: your loop-modelled 3POV from the MD setup
load <YOUR_3POV>.pdb, sox

# AMBER names back to standard
alter resn HIE+HID+HIP, resn='HIS'
alter resn CYX, resn='CYS'
alter resn WAT+TIP3, resn='HOH'
alter resn MG2, resn='MG'
sort

# drop DNA, hydrogens, salt ions, altlocs
remove polymer.nucleic
remove resn DA+DT+DG+DC+DA5+DA3+DT5+DT3+DG5+DG3+DC5+DC3
remove hydro
remove resn NA+CL+K
remove not alt ''+A
alter all, alt=''

# keep only waters coordinating the Mg (first shell)
remove solvent and not byres (solvent within 2.8 of resn MG)

# chain IDs MODELLER can split on: protein A, Mg M, waters W
alter polymer.protein, chain='A'
alter resn MG, chain='M'
alter solvent, chain='W'
sort

print("Mg:", cmd.count_atoms("resn MG"), " waters:", cmd.count_atoms("solvent and elem O"))
save 3pov_tmpl.pdb
