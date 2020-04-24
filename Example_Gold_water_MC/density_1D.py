from MDAnalysis import *
import numpy as np
import matplotlib.pyplot as plt

u = Universe("system_init.pdb", "equil_MC.dcd")

startFrame = 700
frameCount = u.trajectory.n_frames - startFrame
# this is in case of a long trajectory, for this 
# trajectory we can average every couple timesteps
avgfreq=1
nFramesToAverage = int(frameCount / avgfreq )

# these are for calculating electrode separation
index_electrode1=362
index_electrode2=363

bins = 100
h2o  = u.select_atoms("resname HOH and name O")
all_atoms = u.select_atoms("all")

#u.trajectory[0]
#print('z positions of electrodes ' )
#print( all_atoms.positions[index_electrode1][2] , all_atoms.positions[index_electrode2][2] )

counts = [0 for y in range(bins)]

for i in range(nFramesToAverage):
    currentFrame = startFrame + i * avgfreq
    if currentFrame >= u.trajectory.n_frames:
        break
    u.trajectory[currentFrame]

    Lcell = all_atoms.positions[index_electrode2][2] - all_atoms.positions[index_electrode1][2]

    dz = Lcell / bins
    for atom in h2o.positions:
        # position relative to left electrode
        L_relative = atom[2] - all_atoms.positions[index_electrode1][2]
        counts[int(L_relative / dz)] += 1

counts=np.asfarray(counts)

# normalize
counts = counts / len(h2o) / nFramesToAverage / bins

for i in range(len(counts)):
    print( i , counts[i] )

