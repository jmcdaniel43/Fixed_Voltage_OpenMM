from MDAnalysis import *
import numpy as np
from MDAnalysis.core.universe import Merge

u=Universe('nanotube.pdb', 'nanotube.pdb')
atoms=u.select_atoms('resname equals nanc')

external_charge_distance = 20.0

# get center of nanotube
atoms=u.select_atoms('resname equals nanc')
atoms.centroid()
# xmin, xmax of nanotube
x_min = np.amin( atoms.positions[:,0] )
x_max = np.amax( atoms.positions[:,0] )
# radius of nanotube
radius = ( np.amax( atoms.positions[:,1] ) - np.amin( atoms.positions[:,1] ) ) / 2.0

#*********
# creates a new bare-bones universe ...
#*********
def create_new_universe( n_residues , atoms_per_residue , atom_name , positions , resid ):

    n_atoms = n_residues * atoms_per_residue
    resindices = np.repeat(range(n_residues), atoms_per_residue)
    segindices = [0] * n_residues
    u_new = Universe.empty(n_atoms,
                    n_residues=n_residues,
                    atom_resindex=resindices,
                    residue_segindex=segindices,
                    trajectory=True) # necessary for adding coordinates
    u_new.add_TopologyAttr('name', [atom_name,]*n_residues*atoms_per_residue)
    u_new.add_TopologyAttr('type', [atom_name,]*n_residues*atoms_per_residue)
    u_new.add_TopologyAttr('resname', [atom_name, ]*n_residues)
    u_new.add_TopologyAttr('resid', list(range(resid, n_residues+resid)) )
    u_new.atoms.positions = positions

    return u_new


def create_line_of_atoms( center_in , shift , tangent , atom_min , atom_max, number_atoms , atom_name , resid ):
    # get rid of tangent component of center
    center = center_in -  tangent * np.dot( center_in , tangent )
    positions = []
    dr = float( atom_max - atom_min ) / float(number_atoms)
    # now add atoms along specified line
    for i in range( number_atoms ):
        atom_vec = center + shift + tangent * ( atom_min + i * dr )
        positions.append( atom_vec )

    u_new = create_new_universe( 1 , number_atoms , atom_name , positions , resid )
    return u_new

# external line charge
u2 = create_line_of_atoms( atoms.centroid() , np.array( [ 0.0 , 0.0, external_charge_distance ] ) , np.array( [ 1.0 , 0.0 , 0.0 ] ) , x_min , x_max , 40 , 'He' , 2 )

# image line charge at center of cylindar
u3 = create_line_of_atoms( atoms.centroid() , np.array( [ 0.0 , 0.0 , 0.0 ] ) , np.array( [ 1.0 , 0.0 , 0.0 ] ) , x_min , x_max , 40 , 'He' , 3 )

# image line charge at radius^2 / external_charge_distance from center
displacement = radius**2/external_charge_distance
u4 = create_line_of_atoms( atoms.centroid() , np.array( [ 0.0 , 0.0 , displacement ] ) , np.array( [ 1.0 , 0.0 , 0.0 ] ) , x_min , x_max , 40 , 'Ne' , 4 )

atoms2 = u2.select_atoms('all')
atoms3 = u3.select_atoms('all')
atoms4 = u4.select_atoms('all')

# combine nanotube and line charges in new Universe
u5=Merge( atoms , atoms2 , atoms3 , atoms4 )

# now print
atomsprint = u5.select_atoms('all')
atomsprint.write('combine.pdb')
