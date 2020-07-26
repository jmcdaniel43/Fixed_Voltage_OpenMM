from __future__ import print_function
import sys
sys.path.append('../lib/')
#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from MM_classes import *
#*********** Fixed Voltage routines
from Fixed_Voltage_routines import *
#***************************
import numpy as np
import urllib.request
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

# a few run control settings...   WARNING:  write_charge = True will generate a lot of data !!!
simulation_time_ns = 5 ; freq_charge_update_fs = 50 ; freq_traj_output_ps = 5 ; write_charges = True

#******************************************************************
#                Choose type of simulation
simulation_type = "MC_equil"  # either "Constant_V" or "MC_equil"
#**********************************************************************

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# set applied voltage in Volts.  This won't be used for "MC_equil" simulation...
Voltage = 4.0 # in Volts, units will be internally converted later...

# electrode names used to exclude intra-electrode non-bonded interactions ...
#cathode_name="cath"; anode_name = "anod"
# use chain indices instead of residue names to identify electrode...
cathode_index=(2,1,0); anode_index=(3,4,5) # note chain indices start at 0 ...

# Initialize: Input list of pdb and xml files
MMsys=MM( pdb_list = [ 'system_init.pdb', ] , residue_xml_list = [ 'gold_residue.xml' , 'gold-water_residue.xml' ] , ff_xml_list = [ 'gold.xml', 'gold-water.xml' ]  )

# if periodic residue, call this
MMsys.set_periodic_residue(True)

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform('OpenCL')   # only 'Reference' platform is currently implemented!

# initialize Virtual Electrodes, these are electrode `sheets' that solve electrostatics for constant Voltage ...
# can choose electrodes by residue name (this is default)
# can also choose electrodes by chain name (set chain=True)
# can input tuple exclude_element with elements to exclude from virtual electrode, such as dummy Hydrogen atoms ...
MMsys.initialize_electrodes( Voltage, cathode_identifier = cathode_index , anode_identifier = anode_index , chain=True , exclude_element=("H",) )  # chain indices instead of residue names

# initialize atoms indices of electrolyte, we need this for analytic charge correction.  Currently we electrode residue > 100 atoms, electrolyte residue < 100 atoms ... this should be fine?
MMsys.initialize_electrolyte(Natom_cutoff=100)  # make sure all electrode residues have greater than, and all electrolyte residues have less than this number of atoms...

# IMPORTANT: generate exclusions for SAPT-FF.  If flag_SAPT_FF_exclusions = True , then will assume SAPT-FF force field and put in appropriate exclusions.
# set flag_SAPT_FF_exclusions = False if not using SAPT-FF force field
MMsys.generate_exclusions( water_name = 'HOH' , flag_hybrid_water_model = True , flag_SAPT_FF_exclusions = False  )

state = MMsys.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
positions = state.getPositions()

print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))
for j in range(MMsys.system.getNumForces()):
    f = MMsys.system.getForce(j)
    print(type(f), str(MMsys.simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))



# write initial pdb with Drudes, and setup trajectory output
PDBFile.writeFile(MMsys.simmd.topology, positions, open('start_drudes.pdb', 'w'))

if simulation_type == "MC_equil":
    # Monte Carlo equilibration, initialize parameters ...  currently set to move Anode, keep Cathode fixed...
    celldim = MMsys.simmd.topology.getUnitCellDimensions()
    MMsys.MC = MC_parameters(  MMsys.temperature , celldim , electrode_move="Anode" , pressure = 1.0*bar , barofreq = 100 , shiftscale = 0.2 )
    trajectory_file_name = 'equil_MC.dcd'
else :
    trajectory_file_name = 'FV_NVT.dcd'

MMsys.set_trajectory_output( trajectory_file_name , freq_traj_output_ps * 1000 )


for i in range( int(simulation_time_ns * 1000 / freq_traj_output_ps ) ):
    state = MMsys.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))
    for j in range(MMsys.system.getNumForces()):
        f = MMsys.system.getForce(j)
        print(type(f), str(MMsys.simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))


    #**********  Monte Carlo Simulation ********
    if simulation_type == "MC_equil":
        for j in range( int(freq_traj_output_ps * 1000 / MMsys.MC.barofreq) ):
            MMsys.MC_Barostat_step()

    #**********  Constant Voltage Simulation ****
    elif simulation_type == "Constant_V":
        for j in range( int(freq_traj_output_ps * 1000 / freq_charge_update_fs) ):
            # Fixed Voltage Electrostatics ..
            MMsys.Poisson_solver_fixed_voltage( Niterations=4 )
            MMsys.simmd.step( freq_charge_update_fs )
        if write_charges :
            # write charges...
            MMsys.write_electrode_charges( chargeFile )

    else:
        print('simulation type not recognized ...')
        sys.exit()

print('done!')
sys.exit()





