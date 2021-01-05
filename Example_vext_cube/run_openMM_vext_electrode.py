from __future__ import print_function
import sys
sys.path.append('../lib/')
sys.path.append('../../QM_MM/lib/')
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

parser = argparse.ArgumentParser()
parser.add_argument("job", type=str, help="Directory with specific job info")
args = parser.parse_args()

jobname = import_module( args.job )

# a few run control settings...   WARNING:  write_charge = True will generate a lot of data !!!
simulation_time_ns = 5 ; freq_charge_update_fs = 50 ; freq_traj_output_ps = 50 ; write_charges = True

chargeFile = open("charges.dat", "w")
Potential_cubefile = open( "potential.cube" , "w")

#******************************************************************
#                Choose type of simulation
simulation_type = "Constant_V"  # either "Constant_V" or "MC_equil"
#**********************************************************************

ffdir='/home/mcdanielgroup/data/Jesse/Fixed_Voltage_OpenMM/Example_Electrodes_nanotube_BMIM_BF4_ACNT_FV/electrode_ffdir/'

#************************** download SAPT-FF force field files from github
url1 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt.xml"
url2 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt_residues.xml"
filename1, headers = urllib.request.urlretrieve(url1, "sapt.xml")
filename2, headers = urllib.request.urlretrieve(url2, "sapt_residues.xml")

# add extra CustomNonbonded force parameters to .xml file
add_CustomNonbondedForce_SAPTFF_parameters( xml_base = "sapt.xml" , xml_param = ffdir + "graph_customnonbonded.xml" , xml_combine = "sapt_add.xml" )


# *********************************************************************
#                     Create MM system object
#**********************************************************************

Voltage = 2.0
nm_to_bohr = 18.89726

pme_grid_size_a=200
pme_grid_size_b=200
pme_grid_size_c=200
pme_alpha = 7.0

# Initialize: Input list of pdb and xml files, and QMregion_list
MMsys=MM_FixedVoltage( pdb_list = [ jobname.pdbfile, ] , residue_xml_list = [ 'sapt_residues.xml' , ffdir + 'graph_residue_c.xml', ffdir + 'nanotube9x9_residue_c.xml', ffdir + 'graph_residue_n.xml', ffdir + 'nanotube9x9_residue_n.xml' ] , ff_xml_list = [ 'sapt_add.xml', ffdir + 'graph.xml', ffdir + 'graph_c_freeze.xml', ffdir + 'nanotube9x9_c_freeze.xml' , ffdir + 'graph_n_freeze.xml', ffdir + 'nanotube9x9_n_freeze.xml' ]  )
# second MMsys on Reference platform to compute vext
MMsys_vext=MM_FixedVoltage_QMMM( pdb_list = [ jobname.pdbfile, ] , residue_xml_list = [ 'sapt_residues.xml' , ffdir + 'graph_residue_c.xml', ffdir + 'nanotube9x9_residue_c.xml', ffdir + 'graph_residue_n.xml', ffdir + 'nanotube9x9_residue_n.xml' ] , ff_xml_list = [ 'sapt_add.xml', ffdir + 'graph.xml', ffdir + 'graph_c_freeze.xml', ffdir + 'nanotube9x9_c_freeze.xml' , ffdir + 'graph_n_freeze.xml', ffdir + 'nanotube9x9_n_freeze.xml' ] ,  qmmm_ewald = 'True' , pme_alpha= pme_alpha , pme_grid_size = pme_grid_size_a )

MMsys_vext.nbondedForce.setPMEParameters( pme_alpha , pme_grid_size_a , pme_grid_size_b , pme_grid_size_c ) # alpha, nx, ny, nz

# if periodic residue, call this
MMsys.set_periodic_residue(True)
MMsys_vext.set_periodic_residue(True)

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform('OpenCL')   # only 'Reference' platform is currently implemented!
MMsys_vext.set_platform('Reference')


# initialize Virtual Electrodes, these are electrode `sheets' that solve electrostatics for constant Voltage ...
# can choose electrodes by residue name (this is default)
# can also choose electrodes by chain name (set chain=True)
# can input tuple exclude_element with elements to exclude from virtual electrode, such as dummy Hydrogen atoms ...
MMsys.initialize_electrodes( Voltage, cathode_identifier = jobname.nanotube_data.cathode_index , anode_identifier = jobname.nanotube_data.anode_index , chain=True , exclude_element=("H",) , NanoTubes=jobname.nanotube_data.NanoTubes, nanotube_axis=jobname.nanotube_data.nanotube_axis )  # chain indices instead of residue names
MMsys_vext.initialize_electrodes( Voltage, cathode_identifier = jobname.nanotube_data.cathode_index , anode_identifier = jobname.nanotube_data.anode_index , chain=True , exclude_element=("H",) , NanoTubes=jobname.nanotube_data.NanoTubes, nanotube_axis=jobname.nanotube_data.nanotube_axis )  # chain indices instead of residue names


# initialize atoms indices of electrolyte, we need this for analytic charge correction.  Currently we electrode residue > 100 atoms, electrolyte residue < 100 atoms ... this should be fine?
MMsys.initialize_electrolyte(Natom_cutoff=100)  # make sure all electrode residues have greater than, and all electrolyte residues have less than this number of atoms...
MMsys_vext.initialize_electrolyte(Natom_cutoff=100)  # make sure all electrode residues have greater than, and all electrolyte residues have less than this number of atoms...

# IMPORTANT: generate exclusions for SAPT-FF.  If flag_SAPT_FF_exclusions = True , then will assume SAPT-FF force field and put in appropriate exclusions.
# set flag_SAPT_FF_exclusions = False if not using SAPT-FF force field
MMsys.generate_exclusions( flag_SAPT_FF_exclusions = False )
MMsys_vext.generate_exclusions( flag_SAPT_FF_exclusions = False )

# setting compute_intermediate_forces = True will accelerate convergence, but is expensive ...
MMsys.Poisson_solver_fixed_voltage( Niterations=30 , compute_intermediate_forces = True )

state = MMsys.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
positions = state.getPositions()

print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))
for j in range(MMsys.system.getNumForces()):
    f = MMsys.system.getForce(j)
    print(type(f), str(MMsys.simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))


# copy data from GPU context to CPU context ...
print('starting copy from GPU to CPU ...')
copy_positions_charges_between_contexts( MMsys.simmd , MMsys_vext.simmd , MMsys.nbondedForce , MMsys_vext.nbondedForce )

# get vext on PME grid
state = MMsys_vext.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True )

print('retrieving vext grid')
vext_tot = state.getVext_grid()

#** vext_tot from OpenMM is in kJ/mol/e  units, convert to Volts
vext_tot = np.array( vext_tot ) / 96.485

# get box vectors in bohr
box_bohr = get_box_vectors_Bohr( state , nm_to_bohr )

print('printing cube file...')
# print cube file with vext grid
print_cube_file( Potential_cubefile , vext_tot , box_bohr , pme_grid_size_a , pme_grid_size_b , pme_grid_size_c )

sys.exit()


