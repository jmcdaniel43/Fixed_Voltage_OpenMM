from __future__ import print_function
import sys
sys.path.append('../../lib/')
sys.path.append('../../../QM_MM/lib/')
#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from MM_classes_FV import *
from MM_classes_FV_QMMM import *
#*********** Fixed Voltage routines
from Fixed_Voltage_routines import *
#***************************
import numpy as np
import urllib.request
from importlib import import_module
import argparse
from add_customnonbond_xml import add_CustomNonbondedForce_SAPTFF_parameters
# other stuff
from sys import stdout
from time import gmtime, strftime
from datetime import datetime

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit(2000)

# *********************************************************************
#                Fixed-Voltage MD code for electrode simulations
#
#    This code has been generalized to work for flat electrodes (e.g. graphene), and electrodes
#    functionalized with conducting objects such as spherical buckyballs or cylindrical nanotubes
#
#    This includes a module for MC equilibration of the initial system to equilibrate the density.
#**********************************************************************

pdbfile = 'combine.pdb'

chargeFile = open("charges.dat", "w")
Potential_cubefile = open( "potential.cube" , "w")

#******************************************************************
#                Choose type of simulation
simulation_type = "Constant_V"  # either "Constant_V" or "MC_equil"
#**********************************************************************


#************************** download SAPT-FF force field files from github
url1 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt.xml"
url2 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt_residues.xml"
filename1, headers = urllib.request.urlretrieve(url1, "sapt.xml")
filename2, headers = urllib.request.urlretrieve(url2, "sapt_residues.xml")

# add extra CustomNonbonded force parameters to .xml file
add_CustomNonbondedForce_SAPTFF_parameters( xml_base = "sapt.xml" , xml_param = "graph_customnonbonded.xml" , xml_combine = "sapt_add.xml" )


# *********************************************************************
#                     Create MM system object
#**********************************************************************

# set applied voltage in Volts
Voltage = 4.0 # in Volts, units will be internally converted later...
nm_to_bohr = 18.89726

pme_alpha = 7.0
pme_grid_size_a=200
pme_grid_size_b=200
pme_grid_size_c=200

# second MMsys on Reference platform to compute vext
MMsys_vext=MM_FixedVoltage_QMMM( pdb_list = [ pdbfile, ] , residue_xml_list = [ 'sapt_residues.xml' , 'nanotube9x9_residue_c.xml', 'nanotube9x9_residue_n.xml' , 'line_charge_residue.xml' ] , ff_xml_list = [ 'sapt_add.xml', 'graph.xml', 'nanotube9x9_c_freeze.xml' , 'nanotube9x9_n_freeze.xml' , 'line_charge.xml' ] , qmmm_ewald = 'True' , pme_alpha = pme_alpha , pme_grid_size = pme_grid_size_a )

MMsys_vext.nbondedForce.setPMEParameters( pme_alpha, pme_grid_size_a , pme_grid_size_b , pme_grid_size_c ) # alpha, nx, ny, nz

# if periodic residue, call this
MMsys_vext.set_periodic_residue(True)

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys_vext.set_platform('Reference')


# get vext on PME grid
state = MMsys_vext.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True)

print('retrieving vext grid')
vext_tot = state.getVext_grid()
print('retrieving PME_grid_positions')
PME_grid_positions = state.getPME_grid_positions()

#** vext_tot from OpenMM is in kJ/mol/e  units, convert to Volts
vext_tot = np.array( vext_tot ) / 96.485
#** PME_grid_positions from OpenMM is in nanometers, convert to Bohr for input to Psi4
PME_grid_positions = np.array( PME_grid_positions ) * nm_to_bohr

# get box vectors in bohr
box_bohr = get_box_vectors_Bohr( state , nm_to_bohr )

print('printing cube file...')
# print cube file with vext grid
print_cube_file( Potential_cubefile , vext_tot , box_bohr , pme_grid_size_a , pme_grid_size_b , pme_grid_size_c )

sys.exit()


