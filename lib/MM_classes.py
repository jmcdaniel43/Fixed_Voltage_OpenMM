from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#******* Fixed voltage routines
from Fixed_Voltage_routines import *
#******* exclusions routines
from electrode_sapt_exclusions import *
#***********
import random
import numpy
import subprocess

#*************************** README  **************************************
#  This module defines classes that are used in a QM/MM simulation
#  Interfacing Psi4  and OpenMM
#
#  Because these codes use different units and datastructures, conversions
#  between these are done within the QM/MM classes, so that input and output
#  passed to/from the classes is given in units/datastructures of the
#  particular libraries (Psi4, OpenMM) that it comes from/goes to. 
#
#**************************************************************************



#*************************************************
# This MM class is meant to be very general/versatile for a range of simulation types.
# Currently, the three custom types of simulations it allows are:
#  1) QM/MM turned on by inputing "QMregion_list"  -- must use compiled version of customized OpenMM code
#  2) Fixed-Voltage MD for Supercapacitors with a variety of electrode types
#  3) Monte Carlo equilibration of liquid/solid interfaces (e.g. electrode/electrolyte)
#
#**************************************************
class MM(object):
    # input to init is 3 lists , list of pdb files, list of residue xml files, list of force field xml files .
    # **kwargs input is used to override default settings ...
    def __init__(self, pdb_list , residue_xml_list , ff_xml_list , **kwargs  ):
          #*************************************
          #  DEFAULT RUN PARAMETERS, to overide defaults input  as **kwargs ...
          #**************************************
          self.temperature = 300*kelvin
          self.temperature_drude = 1*kelvin
          self.friction = 1/picosecond
          self.friction_drude = 1/picosecond
          self.timestep = 0.001*picoseconds
          self.small_threshold = 1e-6  # threshold for charge magnitude
          self.cutoff = 1.4*nanometer  
          self.qmmm_ewald = False

          # override default settings if input to **kwargs
          if 'temperature' in kwargs :
              self.temperature = kwargs['temperature']
          if 'cutoff' in kwargs :
              self.cutoff = kwargs['cutoff']              
          if 'qmmm_cutoff' in kwargs :
              self.qmmm_cutoff = kwargs['qmmm_cutoff']              
          if 'QMatoms_list' in kwargs :
              self.QMatoms_list = kwargs['QMatoms_list'] 

          # Check if we are doing QM/MM simulation ...
          if 'qmmm_ewald' in kwargs :
              self.qmmm_ewald = kwargs['qmmm_ewald']

          # load bond definitions before creating pdb object (which calls createStandardBonds() internally upon __init__).  Note that loadBondDefinitions is a static method
          # of Topology, so even though PDBFile creates its own topology object, these bond definitions will be applied...
          for residue_file in residue_xml_list:
               Topology().loadBondDefinitions(residue_file)

          # now create pdb object, use first pdb file input
          self.pdb = PDBFile( pdb_list[0] )

          # create modeller
          self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
          # create force field
          self.forcefield = ForceField(*ff_xml_list)
          # add extra particles
          self.modeller.addExtraParticles(self.forcefield)

          # If QM/MM, add QMregion to topology for exclusion in vext calculation...
          #if self.QMMM :
          #    self.modeller.topology.addQMatoms( self.QMregion_list )


          # polarizable simulation?  Figure this out by seeing if we've added any Drude particles ...
          self.polarization = True
          if self.pdb.topology.getNumAtoms() == self.modeller.topology.getNumAtoms():
              self.polarization = False

          if self.polarization :
              #************** Polarizable simulation, use Drude integrator with standard settings
              self.integrator = DrudeLangevinIntegrator(self.temperature, self.friction, self.temperature_drude, self.friction_drude, self.timestep)
              # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
              self.integrator.setMaxDrudeDistance(0.02)
          else :
              #************** Non-polarizable simulation
              self.integrator = LangevinIntegrator(self.temperature, self.friction, self.timestep)


          # create openMM system object
          self.system = self.forcefield.createSystem(self.modeller.topology, nonbondedCutoff=self.cutoff, constraints=HBonds, rigidWater=True)
          # get force types and set method
          self.nbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == NonbondedForce][0]
          
          self.customNonbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomNonbondedForce][0]
          if self.polarization :
              self.drudeForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == DrudeForce][0]
              # will only have this for certain molecules
              self.custombond = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomBondForce][0]

          # set long-range interaction method
          self.nbondedForce.setNonbondedMethod(NonbondedForce.PME)
          self.customNonbondedForce.setNonbondedMethod(min(self.nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))

    def set_QMregion_parameters( self , QMatoms_list , qmmm_cutoff ):
          self.qmmm_cutoff = qmmm_cutoff
          self.QMatoms_list = QMatoms_list

    def set_trajectory_output( self , filename , write_frequency ):
          self.simmd.reporters = []
          self.simmd.reporters.append(DCDReporter(filename, write_frequency))


    # this sets the force groups to used PBC
    def set_periodic_residue(self, flag):
          for i in range(self.system.getNumForces()):
               f = self.system.getForce(i)
               f.setForceGroup(i)
               # if using PBC
               if flag:
                      # Here we are adding periodic boundaries to intra-molecular interactions.  Note that DrudeForce does not have this attribute, and
                      # so if we want to use thole screening for graphite sheets we might have to implement periodic boundaries for this force type
                      if type(f) == HarmonicBondForce or type(f) == HarmonicAngleForce or type(f) == PeriodicTorsionForce or type(f) == RBTorsionForce:
                            f.setUsesPeriodicBoundaryConditions(True)
                            f.usesPeriodicBoundaryConditions()

    # this sets the PME parameters in OpenMM.  The grid size is important for the accuracy of the external potential
    # in the DFT quadrature, since this is interpolated from the PME grid
    def set_PMEParameters( self , pme_alpha , pme_grid_a , pme_grid_b , pme_grid_c ):
        self.nbondedForce.setPMEParameters( pme_alpha , pme_grid_a , pme_grid_b , pme_grid_c )


    # this sets the platform for OpenMM simulation and initializes simulation object
    #*********** Currently can only use 'Reference' for QM/MM ...
    def set_platform( self, platformname ):
          if platformname == 'Reference':
              self.platform = Platform.getPlatformByName('Reference')
              if self.qmmm_ewald :
                  self.properties = {'ReferenceVextGrid': 'true'}
                  self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
              else :
                  self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
          elif platformname == 'CPU':
              self.platform = Platform.getPlatformByName('CPU')
              if self.qmmm_ewald :
                  self.properties = {'ReferenceVextGrid': 'true'}
                  self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
              else :
                   self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
          elif platformname == 'OpenCL':
              self.platform = Platform.getPlatformByName('OpenCL')
              if self.qmmm_ewald :
                  print( 'Can only run QM/MM simulation with reference platform !')
                  sys.exit()
              else :
                  # we found weird bug with 'mixed' precision on OpenCL related to updating parameters in context for gold/water simulation...
                  #self.properties = {'OpenCLPrecision': 'mixed'} 
                  self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
          elif platformname == 'CUDA':
              self.platform = Platform.getPlatformByName('CUDA')
              self.properties = {'Precision': 'mixed'}
              if self.qmmm_ewald :
                  print( 'Can only run QM/MM simulation with reference platform !')
                  sys.exit()
              else :
                  self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
          else:
              print(' Could not recognize platform selection ... ')
              sys.exit(0)
          self.simmd.context.setPositions(self.modeller.positions)



    #***********************************************
    # this initializes the Electrode objects for Constant-Voltage simulation
    # if input chain = True, we initialize electrodes by chain rather than residue name
    def initialize_electrodes( self, Voltage, cathode_identifier , anode_identifier , chain=False, exclude_element=(), **kwargs ):
        # first create electrode objects
        self.Cathode = Electrode_Virtual( cathode_identifier , "cathode" , Voltage , self , chain , exclude_element )
        self.Anode   = Electrode_Virtual( anode_identifier   , "anode"   , Voltage , self , chain , exclude_element )

        # add any Conductors on electrodes, currently, these could be "Buckyballs" or "Nanotubes" ...
        # FIX! Assume all Conductors are on Cathode, see below code.  Need to generalize this!
        self.Conductor_list = []
        if 'BuckyBalls' in kwargs :
            list_temp = kwargs['BuckyBalls'] # this is a list of identifiers (residue or chain) of BuckyBalls
            for identifier in list_temp:
                Buckyball = Buckyball_Virtual( identifier , "cathode" , Voltage , self, chain , exclude_element ) 
                self.Conductor_list.append( Buckyball )

        # for now, need to input cylindrical axis of nanotube.  Eventually, we should write code to
        # determine this automatically....
        if 'NanoTubes' in kwargs :
            list_temp = kwargs['NanoTubes'] # this is a list of identifiers (residue or chain) of NanoTubes
            list_temp2 = kwargs['nanotube_axis'] # corresponding list of nanotube axis
            for identifier , nanotube_axis in zip(list_temp, list_temp2):
                # make sure we have nanotube axis
                if nanotube_axis:
                    Nanotube = Nanotube_Virtual( identifier , "cathode" , Voltage , self, chain , exclude_element, nanotube_axis )
                    self.Conductor_list.append( Nanotube )
                else:
                    print('must input nanotube_axis for all nanotubes!')
                    sys.exit()


        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True) 
        positions = state.getPositions()
        boxVecs = self.simmd.topology.getPeriodicBoxVectors()
        # set electrochemical cell parameters...
        self.set_electrochemical_cell_parameters( positions, boxVecs )

        # now initialize charge on the electrodes based on applied Voltage ...
        self.Cathode.initialize_Charge( self.Lgap, self.Lcell, self )
        self.Anode.initialize_Charge( self.Lgap, self.Lcell, self )


    #******************************************
    # this resets the geometry parameters of the electrochemical cell,
    # specifically 'Lcell' and 'Lgap' which are used in the Poisson solver
    # we make this a standalone method rather than putting it in the 'initialize_electrode'
    # method because we need to call externally if we are doing MC moves on electrodes ...
    #******************************************
    def set_electrochemical_cell_parameters( self , positions , boxVecs ):
       
        # need to figure out Lcell, Lgap .  Get Lcell from first atoms in each electrode, assume electrodes are separated along z-axis
        
        # use z coord of 1st atom in each electrode to compute Lcell
        atom_index_cathode = self.Cathode.electrode_atoms[0].atom_index
        atom_index_anode   = self.Anode.electrode_atoms[0].atom_index

        # set these z coordinates, need these for analytic charge evaluation ...
        z_cath = positions[atom_index_cathode][2] / nanometer
        self.Cathode.set_z_pos(z_cath)
        z_anod = positions[atom_index_anode][2] / nanometer
        self.Anode.set_z_pos(z_anod)
        self.Lcell = abs(z_cath - z_anod)

        # now vacuum gap, = full z length of box minus Lcell
        self.Lgap = boxVecs[2][2] / nanometer - self.Lcell  # in nanometers ...



    #***********************************************
    # this initializes a list of all electrolyte atoms to use for analytic correction of Poisson solver...
    #  
    # rather than passing/hard-coding in a list of all electrolyte residue names, lets do something simple that should be pretty robust
    # if a residue has > Natom_cutoff number of atoms, its an electrode residue
    # if a residue has < Natom_cutoff number of atoms, its an electrolyte residue
    #  a reasonable choice of Natom_cutoff=100, which i don't think will ever lead to a bug...
    def initialize_electrolyte( self , Natom_cutoff=100):
        # make a set of electrolyte residue names, so that we don't have to keep counting atom numbers...
        electrolyte_names=set()
        # initialize list of electrolyte residue objects /atom indices
        self.electrolyte_residues=[]
        self.electrolyte_atom_indices=[]
        for res in self.simmd.topology.residues():
            if res.name in electrolyte_names:
                # add to electrolyte list
                self.electrolyte_residues.append(res)
                for atom in res._atoms:
                    self.electrolyte_atom_indices.append(atom.index)
            else:
                # this is a new residue name, see if its an electrolyte residue
                natoms = 0
                for a in res._atoms:
                    natoms+=1    
                if natoms < Natom_cutoff:
                    # this is an electrolyte residue
                    self.electrolyte_residues.append(res)
                    electrolyte_names.add( res.name )
                    # add to electrolyte list
                    for atom in res._atoms:
                        self.electrolyte_atom_indices.append(atom.index)



    #************************************************
    # This is the Fixed-Voltage Poisson Solver to optimize charges
    # on the electrode subject to applied voltage ...
    #************************************************
    def Poisson_solver_fixed_voltage(self, Niterations=3):
      
        # if QM/MM , make sure we turn off vext_grid calculation to save time with forces... turn back on after converged
        if self.qmmm_ewald :
            platform=self.simmd.context.getPlatform()
            #print(' property value '  , platform.getPropertyValue( self.simmd.context , 'ReferenceVextGrid') )
            platform.setPropertyValue( self.simmd.context , 'ReferenceVextGrid' , "false" )

        #********* Analytic evaluation of total charge on electrodes based on electrolyte coordinates
        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
        positions = state.getPositions()
        # compute charge for both anode/cathode
        self.Cathode.compute_Electrode_charge_analytic( self , positions , self.Conductor_list, z_opposite = self.Anode.z_pos ) # this is correct, input z_pos(Anode) for cathode charge, and vice-versa...
        self.Anode.compute_Electrode_charge_analytic( self , positions , self.Conductor_list, z_opposite = self.Cathode.z_pos )

        #print(" initial charge on Cathode/Anode " , self.Cathode.get_total_charge() , self.Anode.get_total_charge() )

        # make sure these are equal and opposite, but only if no additional conductors ...
        if (abs(self.Cathode.Q_analytic)-abs(self.Anode.Q_analytic)) > self.small_threshold and not self.Conductor_list :
           print( "Analytic charges on Cathode and Anode don't match...something is wrong!" )
           sys.exit(0)

        #*********  Self-consistently solve for electrode charges that obey Fixed-Voltage boundary condition ...
        for i_iter in range(Niterations):

            # need Efield on all electrode atoms, get this from forces on virtual electrode sheets ...
            state = self.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
            forces = state.getForces()

            # print to see convergence, do we want this?
            #for j in range(self.system.getNumForces()):
            #    f = self.system.getForce(j)
            #    print(type(f), str(self.simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))

            # first solve charges on cathode ...
            # note that Cathode.Voltage and Anode.Voltage are Total Voltage drop, i realize this is slightly misleading...
            for atom in self.Cathode.electrode_atoms:
                index = atom.atom_index
                q_i_old = atom.charge
                # Ez , don't divide by zero!
                Ez_external = ( forces[index][2]._value / q_i_old ) if abs(q_i_old) > (0.9*self.small_threshold) else 0.
                # new charge that satisfies fixed Voltage boundary condition...
                # when we switch to atomic units on the right, sigma/2*epsilon0 becomes 4*pi*sigma/2 , since 4*pi*epsilon0=1 in a.u.
                q_i = 2.0 / ( 4.0 * numpy.pi ) * self.Cathode.area_atom * (self.Cathode.Voltage / self.Lgap + Ez_external) * conversion_KjmolNm_Au
                # don't allow charges to get below small_threshold, otherwise can't compute Efield next iteration, and will get stuck at zero forever ...
                if abs(q_i) < self.small_threshold:
                    q_i = self.small_threshold  # Cathode, make positive
                atom.charge = q_i
                self.nbondedForce.setParticleParameters(index, q_i, 1.0 , 0.0)            

            # now charges on anode ...
            for atom in self.Anode.electrode_atoms:
                index = atom.atom_index
                q_i_old = atom.charge
                # Ez , don't divide by zero!
                Ez_external = ( forces[index][2]._value / q_i_old ) if abs(q_i_old) > (0.9*self.small_threshold) else 0.
                # new charge that satisfies fixed Voltage boundary condition...
                # when we switch to atomic units on the right, sigma/2*epsilon0 becomes 4*pi*sigma/2 , since 4*pi*epsilon0=1 in a.u.
                q_i = -2.0 / ( 4.0 * numpy.pi ) * self.Anode.area_atom * (self.Anode.Voltage / self.Lgap + Ez_external) * conversion_KjmolNm_Au
                # don't allow charges to get below small_threshold, otherwise can't compute Efield next iteration, and will get stuck at zero forever ...
                if abs(q_i) < self.small_threshold:
                    q_i = -1.0 * self.small_threshold  # Anode, make negative
                atom.charge = q_i
                self.nbondedForce.setParticleParameters(index, q_i, 1.0 , 0.0)

            # now charges on any conductors that are near electrodes...
            if self.Conductor_list:
                for Conductor in self.Conductor_list:
                    self.Numerical_charge_Conductor( Conductor , forces )

                self.nbondedForce.updateParametersInContext(self.simmd.context)
                # because conductors within cell are "part of electrolyte" as far as analytic charge formula is concerned, need to recomput analytic charges here...
                self.Cathode.compute_Electrode_charge_analytic( self , positions , self.Conductor_list, z_opposite = self.Anode.z_pos ) # this is correct, input z_pos(Anode) for cathode charge, and vice-versa...
                self.Anode.compute_Electrode_charge_analytic( self , positions , self.Conductor_list, z_opposite = self.Cathode.z_pos )

            # Now scale charges to exact Analytic normalization....
            self.Scale_charges_analytic_general()
            # update charges in context ...
            self.nbondedForce.updateParametersInContext(self.simmd.context)

        # this call is just for printing converged charges ...
        self.Scale_charges_analytic_general( print_flag = True )

        # if QM/MM , turn vext back on ...
        if self.qmmm_ewald :
            #print(' property value '  , platform.getPropertyValue( self.simmd.context , 'ReferenceVextGrid') )
            platform.setPropertyValue( self.simmd.context , 'ReferenceVextGrid' , "true" )




    #***************************************
    # this solves conducting boundary condition for Conductor on electrode, given electric field on each atom and applied Voltage
    #
    # so far this works for either Buckyballs (spheres) / Nanotubes (cylindars)
    #
    # 'forces' contain all the forces in the system, used for calculating electric field on virtual atoms
    # these should be updated with current charge distribution on conductors
    #
    # 'Conductor' is a conductor object
    #***************************************
    def Numerical_charge_Conductor( self, Conductor, forces ):
       
        #****************************************************************************
        # Step 1:  Image charges on Conductor.  Project Efield to surface normal vector
        #          solve for the image charge on the Conductor such that the normal field
        #          component is zero inside Conductor
        #******************************************************************************

        # Images charges are set on 'Virtual' atoms of Conductor ...
        for atom in Conductor.electrode_atoms:
            index = atom.atom_index
            (q_i_quantity, sig, eps) = self.nbondedForce.getParticleParameters(index)
            q_i = q_i_quantity._value # quantity = value * units ...

            E_external=[]
            # normal component of Field...
            if abs(q_i) > (0.9*self.small_threshold): 
                E_external.append( forces[index][0]._value / q_i ) # Ex
                E_external.append( forces[index][1]._value / q_i ) # Ey
                E_external.append( forces[index][2]._value / q_i ) # Ez

                # project out normal
                En_external = numpy.dot( numpy.array( E_external ) , numpy.array( [ atom.nx , atom.ny , atom.nz ] ) )
                # now solve for surface charge, requiring Enormal be zero inside conductor...
                q_i = 2.0 / ( 4.0 * numpy.pi ) * Conductor.area_atom * En_external * conversion_KjmolNm_Au

                #print( "normal" , atom.nx , atom.ny , atom.nz , En_external , q_i )

            # don't allow charges to stay below small_threshold, otherwise can't compute Efield next iteration, and will get stuck at zero forever ...
            else: 
                q_i = self.small_threshold  # Cathode, make positive

            atom.charge = q_i
            self.nbondedForce.setParticleParameters(index, atom.charge, sig , eps)


        self.nbondedForce.updateParametersInContext(self.simmd.context)
        state = self.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
        forces = state.getForces()


        #****************************************************************************
        # Step 2:  Charge transfer to Conductor.  Distribute uniformly on atoms.
        #          this is determined from electric field to the right (for cathode) of closest electrode atom being zero
        #          so that the Conductor is at the same Potential as the electrode...
        #******************************************************************************

        # index of close contact atom ...
        conductor_atom = Conductor.Electrode_contact_atom
        conductor_atom_index = conductor_atom.atom_index 
        (q_i_quantity, sig, eps) = self.nbondedForce.getParticleParameters(conductor_atom_index)        
        q_i = q_i_quantity._value # quantity = value * units ...

        # get field normal to surface, most likely this will be in Z-direction, but use general code....
        E_external=[]
        # normal component of Field...
        if abs(q_i) > (0.9*self.small_threshold):
            E_external.append( forces[conductor_atom_index][0]._value / q_i ) # Ex
            E_external.append( forces[conductor_atom_index][1]._value / q_i ) # Ey
            E_external.append( forces[conductor_atom_index][2]._value / q_i ) # Ez

            # project out normal
            En_external = numpy.dot( numpy.array( E_external ) , numpy.array( [ conductor_atom.nx , conductor_atom.ny , conductor_atom.nz ] ) )
        else:
            En_external = 0.0


        # the boundary condition depends on whether the contact is with the Electrode with applied Voltage, or another conductor...
        if Conductor.close_conductor_Electrode :
            # Electrostatics must satisfy on L/R of electrode atom:
            # Left:  -dV/L = - sigma/2eps + Eext + dE_conductor
            # Right:    0  =   sigma/2eps + Eext + dE_conductor
            # therefore, sigma/eps = dV/L
            # and dE_conductor = -( Eext + dV/2L )
            dE_conductor = - ( En_external + self.Cathode.Voltage / self.Lgap / 2.0 ) * conversion_KjmolNm_Au
        else :
            # this is another conductor, no explicit delta_V / L ...
            # there can be no surface charge at this element, because E=0 inside and outside the surface for both boundary conditions
            dE_conductor = - En_external * conversion_KjmolNm_Au


        # Charge depends on geometry of conductor, in general, Q = E * A / 4 *pi where A is area of volume in Gauss' Law integration ...
        if type(Conductor).__name__ == "Buckyball_Virtual" :
            # if buckyball is postive z displacement from cathode, then the field points in negative z for positive charge...
            sign=-1.0
            dQ_conductor =  sign * dE_conductor * Conductor.dr_center_contact**2 

        elif type(Conductor).__name__ == "Nanotube_Virtual" :
            sign=-1.0
            dQ_conductor =  sign * dE_conductor * Conductor.dr_center_contact * Conductor.length / 2.0

        else :
            print( "can't recognize Conductor type in Numerical_charge_Conductor method!" )
            sys.exit()


        #print ( 'dQ_conductor' , dQ_conductor )

        # per atom charge
        dq_atom = dQ_conductor / Conductor.Natoms

        # now ADD this excess charge to Conductor
        for atom in Conductor.electrode_atoms:
            index = atom.atom_index
            (q_i_quantity, sig, eps) = self.nbondedForce.getParticleParameters(index)
            q_i = q_i_quantity._value +  dq_atom
            atom.charge = q_i
            self.nbondedForce.setParticleParameters(index, q_i, sig , eps)



    #***************************************
    # this scales charges to analytic normalization, 
    # which depends on the geometry of the electrodes.
    #
    # for two flat electrodes, the analytic normalization is
    # done independently and this method merely serves as a wrapper
    #
    # however if we have additional conductors (Buckyballs, nanotubes) at
    # the electrodes, then becomes more complicated ...
    #******************************************
    def Scale_charges_analytic_general(self , print_flag = False ):

        #NOTE:  NEED to Generalize.  Currently, assume Conductors are on Cathode if they are present

        if self.Conductor_list:
           # assume anode is scaled normally...
           self.Anode.Scale_charges_analytic( self , print_flag )
           # get analytic correction from anode
           Q_analytic = -1.0 * self.Anode.Q_analytic
            
           # add up total charge on Cathode and all conductors
           Q_numeric_total = self.Cathode.get_total_charge() 

           # loop over Conductors ...
           for Conductor in self.Conductor_list:
               # add charges from 'virtual' atoms
               Q_numeric_total += Conductor.get_total_charge()
           
           if print_flag :
               print( "Q_numeric , Q_analytic charges on Cathode and extra conductors" , Q_numeric_total , Q_analytic )

           # scale factor, make sure not to divide by zero on rare occasions ...
           scale_factor = -1
           if abs(Q_numeric_total) > self.small_threshold:
               scale_factor = Q_analytic / Q_numeric_total

           # now scale all charges on Cathode and conductors
           if scale_factor > 0.0:
               # loop over atoms in Cathode
               for atom in self.Cathode.electrode_atoms:
                   atom.charge = atom.charge * scale_factor
                   self.nbondedForce.setParticleParameters(atom.atom_index, atom.charge, 1.0 , 0.0)
               # loop over Conductors
               for Conductor in self.Conductor_list:
                   for atom in Conductor.electrode_atoms:
                       atom.charge = atom.charge * scale_factor
                       self.nbondedForce.setParticleParameters(atom.atom_index, atom.charge, 1.0 , 0.0)

        else:
            # no extra conductors, scale each electrode to individual Analytic normalization....
            self.Cathode.Scale_charges_analytic( self , print_flag )
            self.Anode.Scale_charges_analytic( self , print_flag )




    #***************************************
    # this generates exclusions for intra-electrode interactions,
    # 
    # if flag_SAPT_FF_exclusions=True, then will also set exclusions for SAPT-FF force field...
    #***************************************
    def generate_exclusions(self, water_name = 'HOH', flag_hybrid_water_model = False ,  flag_SAPT_FF_exclusions = True ):
        # first electrodes, make temporary list of electrode atom indices to pass to exclusions subroutine
        cathode_list=[]
        for atom in self.Cathode.electrode_atoms:
            cathode_list.append( atom.atom_index )
        anode_list=[]
        for atom in self.Anode.electrode_atoms:
            anode_list.append( atom.atom_index )

        # first electrostatic exclusions between all atoms in principle electrode sheet
        exclusion_Electrode_NonbondedForce(self.simmd , self.system, cathode_list, cathode_list, self.customNonbondedForce , self.nbondedForce )
        exclusion_Electrode_NonbondedForce(self.simmd , self.system, anode_list, anode_list, self.customNonbondedForce , self.nbondedForce)

        # now see if we need to add exclusions between any other chains in Cathode
        if len( self.Cathode.electrode_extra_exclusions ) > 0 :
            for chain1 in range(len( self.Cathode.electrode_extra_exclusions )):
                # first exclude between principle electrode sheet and other chains
                exclusion_Electrode_NonbondedForce(self.simmd , self.system, cathode_list, self.Cathode.electrode_extra_exclusions[chain1], self.customNonbondedForce , self.nbondedForce )
                # now between extra chains
                for chain2 in range(chain1 , len( self.Cathode.electrode_extra_exclusions )):
                    exclusion_Electrode_NonbondedForce(self.simmd , self.system, self.Cathode.electrode_extra_exclusions[chain1], self.Cathode.electrode_extra_exclusions[chain2], self.customNonbondedForce , self.nbondedForce )

        # now see if we need to add exclusions between any other chains in Anode
        if len( self.Anode.electrode_extra_exclusions ) > 0 :
            for chain1 in range(len( self.Anode.electrode_extra_exclusions )):
                # first exclude between principle electrode sheet and other chains
                exclusion_Electrode_NonbondedForce(self.simmd , self.system, anode_list, self.Anode.electrode_extra_exclusions[chain1], self.customNonbondedForce , self.nbondedForce )
                # now between extra chains
                for chain2 in range(chain1 , len( self.Anode.electrode_extra_exclusions )):
                    exclusion_Electrode_NonbondedForce(self.simmd , self.system, self.Anode.electrode_extra_exclusions[chain1], self.Anode.electrode_extra_exclusions[chain2], self.customNonbondedForce , self.nbondedForce )

        # now exclusions for Conductors on Electrodes ... DON'T exclude virtual/virtual, exclude real/real and virtual/real
        for Conductor in self.Conductor_list:
            # temporary lists
            Conductor_real_list=[]; Conductor_virtual_list=[]
            for atom in Conductor.electrode_atoms:
                Conductor_virtual_list.append( atom.atom_index )
            for atom in Conductor.electrode_atoms_real:
                Conductor_real_list.append( atom.atom_index )

            exclusion_Electrode_NonbondedForce(self.simmd , self.system, Conductor_real_list, Conductor_real_list, self.customNonbondedForce , self.nbondedForce ) 
            exclusion_Electrode_NonbondedForce(self.simmd , self.system, Conductor_real_list, Conductor_virtual_list, self.customNonbondedForce , self.nbondedForce )


        # if special exclusion for SAPT-FF force field ...
        if flag_SAPT_FF_exclusions:
            SAPT_FF_exclusions( self )

        # if using a hybrid water model, need to create interaction groups for customnonbonded force....
        if flag_hybrid_water_model:
            generate_exclusions_water(self.simmd, self.customNonbondedForce, water_name )

        # having both is redundant, as SAPT-FF already creates interaction groups for water/other
        if flag_SAPT_FF_exclusions and flag_hybrid_water_model:
            print( "redundant settiong of flag_SAPT_FF_exclusions and flag_hybrid_water_model")
            sys.exit()


        # now reinitialize to make sure changes are stored in context
        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
        positions = state.getPositions()
        self.simmd.context.reinitialize()
        self.simmd.context.setPositions(positions)



    

    #*******************************************
    #   this method performs MC/MD steps for moving electrode sheets to
    #   equilibrate the density of the electrolyte
    #
    #         !!! assumes object self.MC exists !!!!
    #
    #   before calling this method, make sure to initialize MC parameters by
    #   construcing self.MC object of class MC_parameters(object):
    #*******************************************
    def MC_Barostat_step( self ):

        # inner functions ...
        def metropolis(pecomp):
            if pecomp < 0.0 * self.MC.RT :
                return True
            elif (random.uniform(0.0,1.0) < numpy.exp(-pecomp/self.MC.RT)):
                return True

        def intra_molecular_vectors( residue_object , pos_ref, positions_array ):
            intra_vec = []
            # loop over atoms in residue
            for atom in residue_object._atoms:
                pos_res_i = positions_array[atom.index]
                vec_i = pos_res_i - pos_ref
                intra_vec.append(numpy.asarray(vec_i))
            return numpy.asarray(intra_vec)


        self.MC.ntrials += 1

        # ************** normal MD steps ****************
        self.simmd.step(self.MC.barofreq)

        # get final positions
        state = self.simmd.context.getState(getEnergy=True, getPositions=True)
        positions = state.getPositions().value_in_unit(nanometer)

        # energy before move
        oldE = state.getPotentialEnergy()
        # store positions
        oldpos = numpy.asarray(positions)
        newpos = numpy.asarray(positions)     

        # now generate trial move
        # randomly choose move distance ... (-1 , 1) Angstrom for now...
        deltalen = self.MC.shiftscale*(random.uniform(0, 1) * 2 - 1)

        # need a reference point for scaling relative positions.  Choose the stationary electrode for this.
        reference_atom_index = -1

        # Currently, can only move Anode, because we might have other conductors on Cathode (Buckyballs, nanotubes)
        # could easily generalize this to move other conductors as well, but haven't yet...
        if self.MC.electrode_move == "Anode" :
            reference_atom_index = self.Cathode.electrode_atoms[0].atom_index # since we are moving Anode, choose stationary cathode as reference...
            # move Anode
            for atom in self.Anode.electrode_atoms:
                newpos[atom.atom_index,2] += deltalen
            # now see if we need to move any other chains in Anode
            if len( self.Anode.electrode_extra_exclusions ) > 0 :
                for extra_anode_sheet in self.Anode.electrode_extra_exclusions :
                    for index in extra_anode_sheet:
                        newpos[index,2] += deltalen
        else:
            print( "Currently, can only move Anode in MC_Barostat_step...need to generalize for other Conductors ...")
            sys.exit()

        Lcell_old = self.Lcell
        Lcell_new = Lcell_old + deltalen


        N_electrolyte_mol=0
        # now loop over electrolyte molecules and move their COM
        for res in self.electrolyte_residues:
            N_electrolyte_mol += 1
            # use first atom in residue as reference...
            for atom in res._atoms:
                pos_ref = newpos[atom.index]
                break
            # get relative coordinates          
            intra_vec = intra_molecular_vectors( res , pos_ref, newpos )

            # convert to coordinate system relative to stationary electrode ...
            pos_ref[2] = pos_ref[2] - newpos[reference_atom_index,2]
            # now scale this reference position by ratio of change in electrochemical cell volume
            pos_ref[2] = pos_ref[2] * Lcell_new / Lcell_old
            # now back to global coordinate system ...
            pos_ref[2] = pos_ref[2] + newpos[reference_atom_index,2]
            # now update positions of all atoms in this molecule
            index_in_molecule = 0
            for atom in res._atoms:
                newpos[atom.index] = pos_ref + intra_vec[index_in_molecule]
                index_in_molecule+=1


        #  Energy of trial move 
        self.simmd.context.setPositions(newpos)
        statenew = self.simmd.context.getState(getEnergy=True,getPositions=True)
        newE = statenew.getPotentialEnergy()
      
        w = newE-oldE + self.MC.pressure*(deltalen * nanometer) - N_electrolyte_mol * self.MC.RT * numpy.log(Lcell_new/Lcell_old)
        if metropolis(w):
            self.MC.naccept += 1
            # move is accepted, update electrochemical cell parameters...
            boxVecs = self.simmd.topology.getPeriodicBoxVectors()
            positions = statenew.getPositions()
            self.set_electrochemical_cell_parameters( positions, boxVecs )
        else:
            # move is rejected, revert to old positions ...
            self.simmd.context.setPositions(oldpos)

        if self.MC.ntrials > 50 :
            print(" After 50 more MC steps ...")
            print("dE, exp(-dE/RT) ", w, numpy.exp(-w/ self.MC.RT))
            print("Accept ratio for last 50 MC moves", self.MC.naccept / self.MC.ntrials)
            if (self.MC.naccept < 0.25*self.MC.ntrials) :
                self.MC.shiftscale /= 1.1
            elif self.MC.naccept > 0.75*self.MC.ntrials :
                self.MC.shiftscale *= 1.1
            # reset ...
            self.MC.ntrials = 0
            self.MC.naccept = 0



    # this method sets an umbrella potential constraining the centroid of input molecule "mol1".
    # this can be done multiple ways, as is controlled by the specific **kwarg passed to the method.
    #      **kwargs:  mol2, atomtype, r0centrold :  constrain to distance from atom "atom" on mol2
    #      **kwargs:  z_global : constrain to absolute z position
    def setumbrella(self, mol1, k , **kwargs ):

        #create mol1 group for centroid
        g1=[]
        for res in self.simmd.topology.residues():
            if res.name == mol1:
                for i in range(len(res._atoms)):
                    g1.append(res._atoms[i].index)
                break

        # option 1: input mol2, atomtype, r0centrold :  constrain to distance from atom "atom" on mol2
        if ('mol2' in kwargs) and ('atomtype' in kwargs) and ('r0centroid' in kwargs ) :
            mol2 = kwargs['mol2'] ; atomtype = kwargs['atomtype'] ; r0centroid = kwargs['r0centroid']
            g2=[]
            for res in self.simmd.topology.residues():
                if res.name == mol2:
                    for i in range(len(res._atoms)):
                        if res._atoms[i].name == atomtype:
                            g2.append(res._atoms[i].index)
                    break

            self.Centroidforce = CustomCentroidBondForce(2,"0.5*k*(distance(g1,g2)-r0centroid)^2")
            self.system.addForce(self.Centroidforce)
            self.Centroidforce.addPerBondParameter("k")
            self.Centroidforce.addPerBondParameter("r0centroid")
            self.Centroidforce.addGroup(g1)
            self.Centroidforce.addGroup(g2)
            bondgroups =[0,1]
            bondparam = [k,r0centroid]
            self.Centroidforce.addBond(bondgroups,bondparam)
            self.Centroidforce.setUsesPeriodicBoundaryConditions(True)
            self.Centroidforce.addGlobalParameter('r0centroid',r0centroid)
            #self.Centroidforce.addEnergyParameterDerivative('r0centroid')

            for i in range(self.system.getNumForces()):
                f = self.system.getForce(i)
                f.setForceGroup(i)

        # option 2: input 'z_global' :  constrain to absolute z distance
        elif 'z_global' in kwargs :
            z_global = kwargs['z_global']
            self.ZForce = CustomExternalForce("0.5 * k * periodicdistance(x,y,z,x,y,z0)^2")
            self.system.addForce(self.ZForce)
            # add particles to force
            for index in g1:
                self.ZForce.addParticle(index)
            self.ZForce.addGlobalParameter('z0', z_global)
            self.ZForce.addGlobalParameter('k', k )

        else:
            print("couldn't recognize **kwargs input in setumbrella method...")
            sys.exit()


        # reinitialize context and set positions
        self.simmd.context.reinitialize()
        self.simmd.context.setPositions(self.modeller.positions)



    #************************************************
    # this method writes electrode charges to output file
    #
    #  FIX:  Not sure the best way to determine order???
    #   we might need to write cathode, conductor , anode charges,
    #   or cathode, anode , conductor charges in either order??
    #   how to automate this??
    #************************************************
    def write_electrode_charges( self, chargeFile ):
        # first cathode then anode charges
        for atom in self.Cathode.electrode_atoms:
            chargeFile.write("{:f} ".format(atom.charge))          
        #for atom in self.Anode.electrode_atoms:
        #    chargeFile.write("{:f} ".format(atom.charge))

        # loop over additional Conductors (buckyballs/nanotubes) if we have them
        for Conductor in self.Conductor_list:
            for atom in Conductor.electrode_atoms:
                chargeFile.write("{:f} ".format(atom.charge))

        for atom in self.Anode.electrode_atoms:
            chargeFile.write("{:f} ".format(atom.charge))


        # write newline to charge file after charge write
        chargeFile.write("\n")
        chargeFile.flush() # flush buffer



    #***************************
    # Getters that are needed if running QM/MM ...
    #***************************


    #**************************
    # input atom_lists is list of lists of atoms
    # returns a list of lists of elements, charges with one-to-one correspondence...
    #**************************
    def get_element_charge_for_atom_lists( self, atom_lists ):

        element_lists=[]
        charge_lists=[]
        # loop over lists in atom_lists , and add list to element_lists , charge_lists
        for atom_list in atom_lists:
            element_list=[]
            charge_list=[]
            # loop over atoms in topology and match atoms from list...
            for atom in self.simmd.topology.atoms():
                # if in atom_list ..
                if atom.index in atom_list:
                    element = atom.element
                    # get atomic charge from force field...
                    (q_i, sig, eps) = self.nbondedForce.getParticleParameters(atom.index)
                    # add to lists
                    if len(dir(element)) == 40:
                         element_list.append( element.symbol )
                    else:
                         element_list.append( atom.name )
                    charge_list.append( q_i._value )

            # now add to element_lists , charge_lists ..
            element_lists.append( element_list )
            charge_lists.append( charge_list )

        return element_lists , charge_lists

  
    #**************************
    # input atom_lists is list of lists of atoms
    # returns a list of lists of positions with one-to-one correspondence...
    #**************************
    def get_positions_for_atom_lists( self , atom_lists ):

        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
        positions = state.getPositions()
        QM_pos = positions[self.QMatoms_list[0]]._value
        box_vectors = [state.getPeriodicBoxVectors()[i]._value for i in range(3)]
        position_lists=[]
        real_position_lists=[]
        # loop over lists in atom_lists , and add list to position_lists
        for atom_list in atom_lists:
            position_list=[]
            real_position_list=[]
            for index in atom_list:
                ind_pos = positions[index]._value
                real_position_list.append( ind_pos )
                r = get_least_mirror_pos( ind_pos , QM_pos , box_vectors)
                position_list.append( [r[i]+QM_pos[i] for i in range(3)] )
            # now add to position_lists ...
            position_lists.append( position_list )
            real_position_lists.append( real_position_list )

        return position_lists , real_position_lists

    #**************************
    # this method gets the list of atoms that are in the QM region, including the specified QM system atoms
    # input is a tuple of atom indices for the QM system and a cutoff distance in nanometers
    #**************************
    def get_QMregion_list(self):
          # converting Angstrom cutoff to nm
          cutoff = self.qmmm_cutoff/10.0
          # getting the atom index for the first atom listed for each residue
          res_atom_ind = []
          res_list = [res for res in self.simmd.topology.residues()]
          for res in res_list:
               res_atom_ind.append(res._atoms[0].index)
          # getting current box size for minimum mirror image calculation
          state = self.simmd.context.getState( getEnergy=False , getForces=False , getVelocities=True , getPositions=True , getParameters=True )
          pos = state.getPositions()
          box_vectors = [state.getPeriodicBoxVectors()[j]._value for j in range(3)]
          QM_pos = [sum( [pos[i]._value[j] for i in self.QMatoms_list] ) / len( self.QMatoms_list ) for j in range(3)]
          # populating QMregion_list and QMdrudes_list
          QMregion_list = []
          QMdrudes_list = []
          for i in range(len(pos)):
               r = get_least_mirror_pos(QM_pos,pos[i]._value,box_vectors)
               dist = sum([r[i]**2 for i in range(3)])**(0.5)
               if dist < cutoff and i in res_atom_ind:
                     ind = res_atom_ind.index(i)
                     if i in self.QMatoms_list:
                         QM_atoms = [atom.index for atom in res_list[ind]._atoms]
                         QM_elements = [atom.element for atom in res_list[ind]._atoms]
                         for j in range(len(QM_atoms)):
                             if hasattr(QM_elements[j],'symbol'):
                                 QMregion_list.append(QM_atoms[j])
                             else:
                                 QMdrudes_list.append(QM_atoms[j])
                     else:
                          QMregion_list.extend([atom.index for atom in res_list[ind]._atoms])
          self.QMregion_list = tuple( QMregion_list )
          self.QMdrudes_list = tuple( QMdrudes_list )

    #********************
    # this method sets exclusions
    #
    #*******************
    def set_QMregion_exclusion(self):
          system = self.simmd.context.getSystem()
          system.clearQMexclude()
          for index in self.QMregion_list:
               system.addQMexclude( index )
          for index in self.QMdrudes_list:
               system.addQMexclude( index )

#****************************************************
# this is a standalone helper method, outside of class
# input is three list vectors
#****************************************************
def get_least_mirror_pos( i_vec, j_vec, box_vec ):
    # getting index of QM atoms in the QMregion_list to avoid self interaction calculation
    r = i_vec - j_vec
    r -= box_vec[2]*math.floor(r[2]/box_vec[2][2]+0.5)
    r -= box_vec[1]*math.floor(r[1]/box_vec[1][1]+0.5)
    r -= box_vec[0]*math.floor(r[0]/box_vec[0][0]+0.5)

    return r

#****************************************************
# this is a small class used for storing Monte Carlo parameters ...
#****************************************************
class MC_parameters(object):
    def __init__( self , temperature , celldim , electrode_move="Anode" , pressure = 1.0*bar , barofreq = 25 , shiftscale = 0.2 ):
        self.RT = BOLTZMANN_CONSTANT_kB * temperature * AVOGADRO_CONSTANT_NA     
        self.pressure = pressure*celldim[0] * celldim[1] * AVOGADRO_CONSTANT_NA # convert pressure to force ...
        self.electrode_move = electrode_move
        self.barofreq = barofreq
        self.shiftscale = shiftscale
        self.ntrials = 0
        self.naccept = 0
