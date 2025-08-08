# DS
# Python script to analyse binding free energy calculations

import os
import numpy as np
import BFEE2.postTreatment as deltaF

def czarfilepath(molname):
	filepath = []
	rmsd = os.path.join(molname, "BFEE/001_RMSDBound/output/abf_1.abf1.czar.pmf")
	filepath.append(rmsd)
	eutheta = os.path.join(molname, "BFEE/002_EulerTheta/output/abf_1.abf1.czar.pmf")
	filepath.append(eutheta)
	euphi = os.path.join(molname, "BFEE/003_EulerPhi/output/abf_1.abf1.czar.pmf")
	filepath.append(euphi)
	eupsi = os.path.join(molname, "BFEE/004_EulerPsi/output/abf_1.abf1.czar.pmf")
	filepath.append(eupsi)
	ptheta = os.path.join(molname, "BFEE/005_PolarTheta/output/abf_1.abf1.czar.pmf")
	filepath.append(ptheta)
	pphi = os.path.join(molname, "BFEE/006_PolarPhi/output/abf_1.abf1.czar.pmf")
	filepath.append(pphi)
	rtrans = os.path.join(molname, "BFEE/007_r/output/formergepmf/merge.abf1.czar.pmf")
	filepath.append(rtrans)
	rmsunbound = os.path.join(molname, "BFEE/008_RMSDUnbound/output/abf_1.abf1.czar.pmf")
	filepath.append(rmsunbound)
	print (filepath)
	return filepath


calculator = deltaF.postTreatment(300., 'namd')
for fnum in range(1,10):
	czarabfpmffiles = czarfilepath("dir-%d" % fnum)
	calculator = deltaF.postTreatment(300., 'namd')
	parameterlist = [10., 0.1, 0.1, 0.1, 0.1, 0.1, 45., 10.] # these are system specific, check origninal paper DOI: 10.1038/s41596-021-00676-1
	result = calculator.geometricBindingFreeEnergy(czarabfpmffiles, parameterlist)
	print(result)
	np.save("deltaG.npy", result)

	
