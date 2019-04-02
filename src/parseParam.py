import numpy as np
import random


#Function to read in reference files
def pairPotRef(pathToPpotsRef, pathToPpotsWrite):
	with open(pathToPpotsRef, "r") as f:  lines = f.readlines()

	spec1 = [] #species, i.e. C, H, etc.
	spec2 = []
	pairPotParam = []

	header1 = lines[0]
	header2 = lines[1]

	holdTmp = lines[0].split()
	numPpots = int(holdTmp[1])

	for i in range(2, numPpots+2):
		holdTmp = lines[i].split()
		spec1.append(holdTmp[0])
		spec2.append(holdTmp[1])
		for j in range(2,12): #12 is total # of param to print out, can be added to function input too
			pairPotParam.append(float(holdTmp[j]))

	#Write out reference parameters to the ppots.nonortho file
	with open(pathToPpotsWrite, "w") as f:
		f.write(header1)
		f.write(header2)

		#Also hold the initialized parameters that are to be changed
		holdInitialPpots = []

		for i in range(0, numPpots):
			f.write("%s %s %f %f %f %f %f %f %f %f %f %f \n" %
						( spec1[i], spec2[i], pairPotParam[0 + 10*(i)], pairPotParam[1 + 10*(i)],
						pairPotParam[2 + 10*(i)], pairPotParam[3 + 10*(i)], pairPotParam[4 + 10*(i)],
						pairPotParam[5 + 10*(i)], pairPotParam[6 + 10*(i)], pairPotParam[7 + 10*(i)],
						pairPotParam[8 + 10*(i)], pairPotParam[9 + 10*(i)]))
			holdInitialPpots.extend([pairPotParam[0 + 10*(i)], pairPotParam[1 + 10*(i)], pairPotParam[2 + 10*(i)], pairPotParam[3 + 10*(i)], pairPotParam[4 + 10*(i)]])

	return holdInitialPpots	#	return header1, header2, numPpots, spec1, spec2, pairPotParam

#Function to read in bond integrals
def bondIntsRef(pathToBondIntsRef, pathToBondIntsWrite):
    with open(pathToBondIntsRef, "r") as f: lines = f.readlines()

    bispec1 = []
    bispec2 = []
    btype = []
    biParam = []
    overlapParam = []
    biheader1 = lines[0]
    biheader2 = lines[1]

    mytmp = lines[0].split()
    noBis = int(mytmp[1])

    for i in range(2, noBis + 2):
		mytmp = lines[i].split()
		bispec1.append(mytmp[0])
		bispec2.append(mytmp[1])
		btype.append(mytmp[2])
		for j in range(3, 11):
			biParam.append(float(mytmp[j])) #bond integrals
		for j in range(11, 19):
			overlapParam.append(float(mytmp[j])) #overlap matrix parameters

	#Write out reference parameters to bondints.nonortho file
    with open(pathToBondIntsWrite, "w") as f:
	    f.write(biheader1)
	    f.write(biheader2)

	    holdInitialBondInts = []
	    for i in range(0, noBis):
			f.write("%s %s %s  %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n" %
						( bispec1[i], bispec2[i], btype[i], biParam[0 + 8*(i)], biParam[1 + 8*(i)],
						biParam[2 + 8*(i)], biParam[3 + 8*(i)], biParam[4 + 8*(i)], biParam[5 + 8*(i)],
						biParam[6 + 8*(i)], biParam[7 + 8*(i)],
						overlapParam[0 + 8*(i)], overlapParam[1 + 8*(i)], overlapParam[2 + 8*(i)],
						overlapParam[3 + 8*(i)], overlapParam[4 + 8*(i)], overlapParam[5 + 8*(i)],
						overlapParam[6 + 8*(i)], overlapParam[7 + 8*(i)]))

			holdInitialBondInts.extend([biParam[0 + 8*(i)], biParam[1 + 8*(i)], biParam[2 + 8*(i)]])


    return holdInitialBondInts
#	return bispec1, bispec2, btype, noBis, biParam, overlapParam

def electronRef(pathToElecRef, pathToElecWrite):
	with open(pathToElecRef, "r") as f: lines = f.readlines()

	eleSpec = []
	eleValence = []
	eleParam = []

	elehead1 = lines[0]
	elehead2 = lines[1]

	mytmp = lines[0].split()
	noEle = int(mytmp[1])

	for i in range(2, noEle + 2):
		mytmp = lines[i].split()
		eleSpec.append(mytmp[0])
		eleValence.append(mytmp[1])
		for j in range(2, 13):
			eleParam.append(float(mytmp[j]))

	#Write out reference parameters to electrons.dat file
	with open(pathToElecWrite, "w") as f:
		f.write(elehead1)
		f.write(elehead2)

		holdInitialElec = []
		for i in range(0, noEle):
			f.write("%s %s %f %f %f %f %f %f %f %f %f %f %f \n" %
						( eleSpec[i], eleValence[i], eleParam[0 + 11*i], eleParam[1 + 11*i],
						 eleParam[2 + 11*i], eleParam[3 + 11*i], eleParam[4 + 11*i], eleParam[5 + 11*i],
						 eleParam[6 + 11*i], eleParam[7 + 11*i], eleParam[8 + 11*i], eleParam[9 + 11*i],
						 eleParam[10 + 11*i]))


			holdInitialElec.extend([eleParam[1 + 11*i], eleParam[2 + 11*i],  eleParam[6 + 11*i]])
	return holdInitialElec
#	return elehead1, elehead2, noEle, eleSpec, eleValence, eleParam

#Concatenate ppots, bondints, elec together
def initializeT(ppotsInitial, bondintsInitial, elecInitial):
	ppots = pairPotRef("ppots.ref", "ppots.nonortho")
	bondints = bondIntsRef("bondints.ref", "bondints.nonortho")
	elec = electronRef("electrons.ref", "electrons.dat")

	initialT = ppots + bondints + elec
	return np.array(initialT)
