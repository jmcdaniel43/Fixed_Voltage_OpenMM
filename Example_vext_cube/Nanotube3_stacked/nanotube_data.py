import numpy as np

# use chain indices instead of residue names to identify electrode...
cathode_index=(0,5); anode_index=(4,9) # note chain indices start at 0 ...
# list of Nanotubes ..
NanoTubes = [(1,6),(2,7),(3,8)]
# for now , must input nanotube axis.  Eventually, write code to automatically determine this...
nanotube_axis=[(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (-0.5, np.sqrt(3)/2.0, 0.0) ]  # axis in cartesian coordinates, not box vectors ...

