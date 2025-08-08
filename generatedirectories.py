# DS
# Python script to generate input for binding free energy calculations

import BFEE2.inputGenerator as igen
import subprocess
import os
from Bio import SeqIO, AlignIO
import numpy as np
import glob

def mymkdir(s):
	if not os.path.exists(s):
		os.makedirs(s)
def mysymlink(source, dest):
	#print source, dest
	if not os.path.exists(dest):
		os.symlink(source, dest)

# this alignment block is optional; see below where the alignment is used
# for selecting specific residues in the ligand.
out_file = "../../../alignment.aln"
align = AlignIO.read(out_file, "clustal")
seqnumber = (str(align[0].seq))
array = np.array([*seqnumber])
newarray=np.nonzero(array != "-")[0][:-1]
alignmentlookup = dict()
for i in range(len(align)):
	print(align[i].seq, align[i].id)
	seqnumber = (str(align[i].seq))
	array = np.array([*seqnumber])
	resids = np.cumsum(array != "-")
	print(resids)
	print(resids[newarray])
	alignmentlookup[int((align[i].id)[3])] = str(resids[newarray])[7:-1].replace(",", "").replace("\n", "").strip()
	#help(align[i])
print(alignmentlookup)
# exit()

# CHARMM36
parameterlist = ["toppar/par_all36m_prot.prm"  ,
"toppar/par_all36_na.prm" ,                     
"toppar/par_all36_lipid.prm" ,                    
"toppar/par_all36_cgenff.prm" ,                 
"toppar/par_RCpigment_modified_heme.inp"  ,                    
"toppar/par_all36_carb.prm"  ,                  
"toppar/par_misc.inp"   ,                  
"toppar/par_all36_lipid_sphingo.prm"   ,                  
"toppar/par_all36_carb_glycolipid.str"  ,                 
"toppar/par_all36_carb_imlab.str"   ,                                      
"toppar/par_water_ions.str"  ,
"toppar/par_ironsulfur.prm"    ]

generator = igen.inputGenerator()
for fnum in range(1,10):
	#subprocess.call("rm -rf dir-%d" % fnum, shell=True)
	mymkdir("dir-%d" % fnum)
	generator.generateNAMDGeometricFiles("dir-%d" % fnum, "dir-%d/system.psf" % fnum, "dir-%d/system.pdb" %fnum, "charmm", parameterlist, 300, "segid A B C", "segid L and resid %s" % alignmentlookup[fnum], membraneProtein=False,stratification=[1, 1, 1, 1, 1, 1, 1, 1], parallelRuns=2)
	
