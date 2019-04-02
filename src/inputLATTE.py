import subprocess
import numpy as np
from updateParam import updatePairPot,updateBondInts,updateElectron


def getNums(fixedParams):
    pathRef = fixedParams[0]
    def countLines(ext):
        with open(pathRef+'/'+ext,'r') as f: return len([l for l in f.readlines() if len(l)>0])-2
    return map(countLines,['ppots.ref','bondints.ref','electrons.ref'])

def tightBind(x,t):
  numPpots,numBondInts,numElec = getNums(x)

  ppots,bondints,elecs = tToTuple(t,numPpots,numBondInts,numElec)

  updatePairPot(ppots)
  updateBondInts(bondints)
  updateElectron(elecs)

  return readDataWriteInput(x) #readDataWriteInput(x)


def tToTuple(t,numPpots,numBondInts,numElec):
    pPots,bondInts,elecs=[],[],[]

    for i in range(numPpots):    pPots.append(t[5*i:(5*i+5)])
    shift = numPpots*5

    for i in range(numBondInts): bondInts.append(t[(shift+3*i):(shift+3*i+3)])
    shift+=numBondInts*3

    for i in range(numElec):     elecs.append(t[(shift+3*i):(shift+3*i+3)])

    return (pPots,bondInts,elecs)




def readDataWriteInput(x): #yangtab, inputblock: paths to these
    nentries   = x[4]
    yangtab    = x[1]
    inputblock = x[2]


    dataFile = open(yangtab, "r")
    lines = dataFile.readlines()
    dataFile.close()

    nats = []
    box = []
    element = []
    coord = []
    force = []
    energy = []
    dipole = []

    count = 0

    for i in range(0, nentries): #loop over all entries
        nats.append(int(lines[count])) # store the number of atoms
        numat = int(lines[count])
        count = count+1

        for j in range(0,3):
            mytmp = lines[count].split()
            for k in range(0,3):
                box.append(float(mytmp[k])) # and the vectors for the pbc
            count = count + 1

        for j in range(0, numat):
            mytmp = lines[count].split()
            count = count + 1
            element.append(mytmp[0]) # Now read the species
            for k in range(1,4):
                coord.append(float(mytmp[k])) # and the coordinates

        energy.append(float(lines[count])) # and the atomization energy
        count = count + 1
        for j in range(0,numat):
            mytmp = lines[count].split()
            count = count + 1
            for k in range(0,3) :
                force.append(float(mytmp[k])) # and finally the forces
        dipole.append(float(lines[count]))
        count = count + 1

    #Loop over entries in yangtab... call function that changes a parameter randomly/decide # of times to loop around?

    refEnergyPerAtom = []
    for i in range(0, nentries):
        refEnergyPerAtom.append(energy[i]/nats[i])

    energyStd = np.std(refEnergyPerAtom)

    #"force" contains all the reference forces (in order: x1, y1, z1, x2, y2, z2, etc.)
    forceStd = np.std(force)

    count = 0

    tightBindingEnergyPerAtom = []
    tightBindingForces = []
    for yangloop in range(0, nentries):
        fcount = count
        inputLATTE = open(inputblock, "w")
        inputLATTE.write("%d\n" % nats[yangloop])
        inputLATTE.write("%14.9f %14.9f %14.9f\n" % (box[0 + yangloop*9], box[1 + yangloop*9], box[2 + yangloop*9]))
        inputLATTE.write("%14.9f %14.9f %14.9f\n" % (box[3 + yangloop*9], box[4 + yangloop*9], box[5 + yangloop*9]))
        inputLATTE.write("%14.9f %14.9f %14.9f\n" % (box[6 + yangloop*9], box[7 + yangloop*9], box[8 + yangloop*9]))

        #Count number of C/H/N/O to subtract off from atomization energy
        numberCarbon = 0
        numberHydrogen = 0
        numberOxygen = 0
        numberNitrogen = 0

        for i in range(0,nats[yangloop]):
            inputLATTE.write("%s %14.9f %14.9f %14.9f \n" % ((element[count], coord[0 + count*3], coord[1 + count*3], coord[2 + count*3])))

            if element[count] == "C": numberCarbon += 1
            if element[count] == "H": numberHydrogen += 1
            if element[count] == "N": numberNitrogen += 1
            if element[count] == "O": numberOxygen += 1

            count = count + 1
        inputLATTE.close()

        with open("lat_tmp.dat","w") as fstore_latte:
            subprocess.call(["./LATTE_DOUBLE"], stdout = fstore_latte)

        latteOutput = open("fittingoutput.dat", "r")
        lines = latteOutput.readlines()
        freeEnergy = float(lines[0])

        for i in range(1, nats[yangloop] + 1):
            mytmp = lines[i].split()
            for j in range(0, 3):
                tightBindingForces.append(float(mytmp[j]))
        #	tbDipole = float(lines[-1])
        latteOutput.close()

        #Energy with isolated atoms subtracted off
        carbonEnergy = -1.2362
        hydrogenEnergy = -1.1170
        nitrogenEnergy = -3.1204
        oxygenEnergy = -1.5153

        freeEnergy = freeEnergy - (float(numberCarbon)*carbonEnergy + float(numberHydrogen)*hydrogenEnergy + float(numberNitrogen)*nitrogenEnergy
                + float(numberOxygen)*oxygenEnergy)

        tightBindingEnergyPerAtom.append(freeEnergy/nats[yangloop])

    return np.array(refEnergyPerAtom), energyStd, np.array(tightBindingEnergyPerAtom), np.array(force), forceStd, np.array(tightBindingForces)
