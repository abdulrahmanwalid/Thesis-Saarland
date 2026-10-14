# Homology models of herpesvirus alkaline nucleases on the KSHV SOX template (3POV)
# Template: 3pov_tmpl.pdb from clean_template.pml (protein + 2 Mg2+ + coordinating waters)
# Fill in the <...> fields below (UniProt IDs are listed in README.md), run once per homolog.

# %% Settings (fill these in)
import os, urllib.request
from modeller import Environ, Model, Alignment
from modeller.automodel import AutoModel, assess

UNIPROT   = "<UNIPROT_ID>"         # target accession, see README.md
NAME      = "<NAME>"               # output name, e.g. UL12, UL98, BGLF5
TMPL_DIR  = "<TEMPLATE_FOLDER>"    # folder containing 3pov_tmpl.pdb
TMPL_CODE = "3pov_tmpl"            # template file name without .pdb
N_MODELS  = 20                     # number of models to build

os.chdir(TMPL_DIR)                 # all MODELLER output lands next to the template

# %% Count Mg2+ and waters in the template (needed for the alignment block)
lines = [l for l in open(TMPL_CODE + ".pdb") if l.startswith(("ATOM", "HETATM"))]
N_MG = sum(1 for l in lines if l[17:20].strip() == "MG")
N_W  = sum(1 for l in lines if l[17:20].strip() == "HOH" and l[12:16].strip() in ("O", "OW"))
print(N_MG, "Mg |", N_W, "waters")

# %% Target sequence from UniProt
fa  = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.fasta").read().decode()
seq = "".join(fa.split("\n")[1:])
print(NAME, len(seq), "aa")        # compare with the length in README.md

# %% Align target to template (Mg + waters carried as block residues)
env = Environ()
env.io.hetatm = True
env.io.water  = True
env.io.atom_files_directory = ["."]

aln = Alignment(env)
aln.append_model(Model(env, file=TMPL_CODE + ".pdb"), align_codes=TMPL_CODE, atom_files=TMPL_CODE)
aln.append_sequence(seq + "/" + "." * N_MG + "/" + "w" * N_W)
aln[1].code = NAME
aln.align2d()
aln.write(file=f"{NAME}.ali", alignment_format="PIR")
aln.write(file=f"{NAME}.pap", alignment_format="PAP")   # check the catalytic residues line up here

# %% Build models and rank by DOPE
a = AutoModel(env, alnfile=f"{NAME}.ali", knowns=TMPL_CODE, sequence=NAME,
              assess_methods=(assess.DOPE,))
a.starting_model, a.ending_model = 1, N_MODELS
a.make()

# DOPE ranking is a starting point only: the final model is chosen by visual inspection (see README.md)
ok = sorted([x for x in a.outputs if x["failure"] is None], key=lambda x: x["DOPE score"])
for x in ok:
    print(x["name"], round(x["DOPE score"], 1))
