#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import *

#*************************** README  **************************************
#  This module defines methods that introduce additional exclusions
#  for electrodes, and specific electrolyte molecules within the SAPT force field.
#  This module should thus be considered "part of" the force field,
#  as (non-electrode) exclusions are force-field specific
#**************************************************************************

"""
    Set up exclusions for electrodes and certain molecules in SAPT force field.
    Exclusions should be implemented separately for different molecule types.
    This SAPT_FF_exclusion class should figure out which methods to subsequently call
    based on the system to generate molecule-specific exclusions
"""


#**********************************
# Adds exclusions for intra-electrode interactions of atoms
# Call this method with lists of electrode atoms 'electrode1' and 'electrode2'
# and all long range interactions between atoms on this electrode(s) will
# be excluded.  These lists could be the same or different, and we take two lists to allow versatility
#**********************************
def exclusion_Electrode_NonbondedForce(sim, system, electrode1, electrode2, customNonbondedForce , nbondedForce ):
    # first figure out which exclusions we already have (1-4 atoms and closer).  The code doesn't
    # allow the same exclusion to be added twice
    flagexclusions = {}

    for i in range(customNonbondedForce.getNumExclusions()):
        (particle1, particle2) = customNonbondedForce.getExclusionParticles(i)
        string1=str(particle1)+"_"+str(particle2)
        string2=str(particle2)+"_"+str(particle1)
        flagexclusions[string1]=1
        flagexclusions[string2]=1

    # now add exclusions for every atom pair in electrode if we don't already have them
    if electrode1 == electrode2:
        for i in range(len(electrode1)):
            indexi = electrode1[i]
            for j in range(i+1,len(electrode2)):
                indexj = electrode2[j]
                string1=str(indexi)+"_"+str(indexj)
                string2=str(indexj)+"_"+str(indexi)
                if string1 in flagexclusions and string2 in flagexclusions:
                    continue
                else:
                    customNonbondedForce.addExclusion(indexi,indexj)
                    nbondedForce.addException(indexi,indexj,0,1,0,True)
    else:
        # different chains/residues
        for i in range(len(electrode1)):
            indexi = electrode1[i]
            for j in range(len(electrode2)):
                indexj = electrode2[j]
                string1=str(indexi)+"_"+str(indexj)
                string2=str(indexj)+"_"+str(indexi)
                if string1 in flagexclusions and string2 in flagexclusions:
                    continue
                else:
                    customNonbondedForce.addExclusion(indexi,indexj)
                    nbondedForce.addException(indexi,indexj,0,1,0,True)





#*********************************
# water models generally use specific water-water interaction parameters
# and separate interaction parameters for water-other
# here assume water-water interactions are in NonbondedForce, and
# water-other are in CustomNonbonded force, and create interaction
# groups for CustomNonbonded force...
#**********************************
def generate_exclusions_water(sim,customNonbondedForce,watername):

    # create Interaction Groups for hybrid water model
    water=set()
    notwater=set()
    print('Creating Interaction Groups for CustomNonBonded.  These interactions will be computed between water-other, not water-water.')
    for res in sim.topology.residues():
        if res.name == watername:
            for i in range(len(res._atoms)):
                water.update([res._atoms[i].index])
        else:
            for i in range(len(res._atoms)):
                notwater.update([res._atoms[i].index])

    customNonbondedForce.addInteractionGroup(water, notwater)
    customNonbondedForce.addInteractionGroup(notwater, notwater)




class SAPT_FF_exclusions(object):
    def __init__(self, MMsys ):

        #************************** List of molecule types that have special exclusions ******
        # now see what molecule types are present.  These name are hardcoded in, i don't see any way around this since special exclusions
        # have to be molecule specific
        self.watername = 'HOH'
        self.TFSIname  = 'Tf2N'
        #***********************************************

        # local copies of force classes
        self.nbondedForce = MMsys.nbondedForce
        self.customNonbondedForce = MMsys.customNonbondedForce 
        self.drudeForce = MMsys.drudeForce
        # will only have this for certain molecules
        self.custombond = MMsys.custombond
  
        # ************* Add Water exclusions
        for res in MMsys.simmd.topology.residues():
            if res.name == self.watername:
                generate_exclusions_water(MMsys.simmd, self.customNonbondedForce, self.watername)
                break

        # ************* Add TFSI exclusions, note might refer to this as Tf2N
        for res in MMsys.simmd.topology.residues():
            if res.name == self.TFSIname:
                self.generate_exclusions_TFSI(MMsys.simmd,MMsys.system)
                break



    def generate_exclusions_TFSI(self,sim,system):
        """
        This creates exclusions for TFSI nonbonded interactions, and update
        Screened Drude interactions.  1-5 non-Coulomb interaction are accounted for
        using CustomBondForce
        """
        print('Creating Exclusions for TFSI')

        # map from global particle index to drudeforce object index
        particleMap = {}
        for i in range(self.drudeForce.getNumParticles()):
            particleMap[self.drudeForce.getParticleParameters(i)[0]] = i

        # can't add duplicate ScreenedPairs, so store what we already have
        flagexceptions = {}
        for i in range(self.nbondedForce.getNumExceptions()):
            (particle1, particle2, charge, sigma, epsilon) = self.nbondedForce.getExceptionParameters(i)
            string1=str(particle1)+"_"+str(particle2)
            string2=str(particle2)+"_"+str(particle1)
            flagexceptions[string1]=1
            flagexceptions[string2]=1

        # can't add duplicate customNonbonded exclusions, so store what we already have
        flagexclusions = {}
        for i in range(self.customNonbondedForce.getNumExclusions()):
            (particle1, particle2) = self.customNonbondedForce.getExclusionParticles(i)
            string1=str(particle1)+"_"+str(particle2)
            string2=str(particle2)+"_"+str(particle1)
            flagexclusions[string1]=1
            flagexclusions[string2]=1

        # add exclusions for all atom pairs on TFSI residues, and when a drude pair is
        # excluded add a corresponding screened thole interaction in its place
        for res in sim.topology.residues():
            if res.name == self.TFSIname:
                for i in range(len(res._atoms)-1):
                    for j in range(i+1,len(res._atoms)):
                        (indi,indj) = (res._atoms[i].index, res._atoms[j].index)
                        # here it doesn't matter if we already have this, since we pass the "True" flag
                        self.nbondedForce.addException(indi,indj,0,1,0,True)
                        # make sure we don't already exlude this customnonbond
                        string1=str(indi)+"_"+str(indj)
                        string2=str(indj)+"_"+str(indi)
                        if string1 in flagexclusions and string2 in flagexclusions:
                            continue
                        else:
                            self.customNonbondedForce.addExclusion(indi,indj)
                        # add thole if we're excluding two drudes
                        if indi in particleMap and indj in particleMap:
                            # make sure we don't already have this screened pair
                            if string1 in flagexceptions or string2 in flagexceptions:
                                continue
                            else:
                                drudei = particleMap[indi]
                                drudej = particleMap[indj]
                                self.drudeForce.addScreenedPair(drudei, drudej, 2.0)




