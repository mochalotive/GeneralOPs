import sys
if len(sys.argv) < 3:
  print("supply (in this order: plumedfile containing -k(phi-phi_t), number of cells, and the resnames in the order in the plumed inputs")
  sys.exit()

forces = []
for resname in range(len(sys.argv)-3):
  forces.append([])
f = open(sys.argv[1])
fl = f.readlines()
fls = []
for i in fl:
  if i[0:6] == "PLUMED":
    continue
  else:
    fls.append(i.strip())
finalForce = fls[(len(sys.argv)-3)*2:len(fls)]
#check for #
for i in range(len(sys.argv)-3):
  if fls[i*2][0:2] != "#!":
    print("malformed plumed output file. Are you sure you're printed all the resname forces?")
    sys.exit()
  else:
    forces[i].append(fls[2*i+1])
numRes = len(sys.argv)- 3
for i in range(len(sys.argv)-3):
  for j in range(len(finalForce)):
    if (j%numRes == i):
      forces[i].append(finalForce[j])
for i in range(len(forces)):
  for j in range(len(forces[i])):
    forces[i][j] = forces[i][j].split()
    forces[i][j] = forces[i][j][1:len(forces[i][j])]
    for k in range(len(forces[i][j])):
      forces[i][j][k] = float(forces[i][j][k])
#print(len(forces))
#print(len(forces[0]))
#print(len(forces[0][0]))
#forces are the gradients of the free energy:
#shape is res X numsteps X numops
#read in op grads from files

#graPhi is gradients of ops wrt atoms
#shape is res X cell X step X atoms
numCell = int(sys.argv[2])

resnames = sys.argv[3:len(sys.argv)]

gradPhi = []
numAtoms = 0
for res in resnames: 
  gradPhi.append([])#res
  for cell in range(numCell):
    gradPhi[-1].append([])#op
    g = open("forces" + res + str(cell))
    gl = g.readlines()
    numAtoms = int(gl[0].strip())
#    print(numAtoms)
    gls = [line.split() for line in gl]
    for atom in range(len(gls)):
      if atom%(numAtoms+2) ==0:
        gradPhi[-1][-1].append([])#step
        continue
      if atom%(numAtoms+2) == 1:
        continue
      else:
        gradPhi[-1][-1][-1].append([float(gls[atom][1]), float(gls[atom][2]), float(gls[atom][3])]) 
#print(len(gradPhi))
#print(len(gradPhi[0]))
#print(len(gradPhi[0][0]))
#print(len(gradPhi[0][0][0]))
#now process data to generate forces w.r.t atoms and make opd file. 

for res in range(len(gradPhi)):
  for cell in range(len(gradPhi[res])):
    for step in range(len(gradPhi[res][cell])):
      for atom in range(len(gradPhi[res][cell][step])):
        if (forces[res][step][cell]) == 0: 
          continue
        gradPhi[res][cell][step][atom][0]/= forces[res][step][cell]
        gradPhi[res][cell][step][atom][1]/= forces[res][step][cell]
        gradPhi[res][cell][step][atom][2]/= forces[res][step][cell]

#write some files out:       
#just plop it all in one file!

jakeFile = open("gradF.csv", "w")
resNum = len(forces)
jakeFile.write(str(len(forces[0])) + "," + str(len(forces)) + "," + str(len(forces[0][0])) + ",")
if (len(gradPhi) != resNum):
  print("inconsistent residue numbers")
  sys.exit()
print(forces)
for step in range(len(forces[0])):
  resNum = len(forces)
  for res in range(resNum):
    cellNum = len(forces[res][step])
    for cell in range(cellNum):
      jakeFile.write(str(forces[res][step][cell]) + "," )
jakeFile.close()
mij = open("matrix.csv", "w")
mij.write(str(len(gradPhi[0][0])) + "," + str(len(gradPhi)) + "," + str(len(gradPhi[0])) + "," + str(len(gradPhi[0][0][0])) + ","   )
for step in range(len(gradPhi[0][0])):
  for res in range(len(gradPhi)):   
    for cell in range(len(gradPhi[res])):
      for atom in range(len(gradPhi[res][cell][step])):
        mij.write(f"{gradPhi[res][cell][step][atom][0]} {gradPhi[res][cell][step][atom][0]} {gradPhi[res][cell][step][atom][0]},")

mij.close()
