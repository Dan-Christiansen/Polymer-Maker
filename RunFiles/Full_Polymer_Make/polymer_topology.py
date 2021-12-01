# Tool to make linear polymer structure and topology files for use with GROMACS OPLS-AA molecular dynamics simulation
#
#
#
# Input files: polymer_formula
# Input arguments:
# Number of polymers in the system
# Rotation angle between neighboring monomers
#
import os
import numpy as np

def oss(fun): # Abbreviated 'os.system' function
    os.system(fun)

def polymer_topology(n_monomers):
    monomer = 'newmol' # Polymer filename, usually newmol

    # Produce a
    os.system('/usr/local/gromacs_gmx/bin/gmx pdb2gmx -f '+monomer+'.gro -ff oplsaa -p topol_mon.top -water spce')

    # Flag for when to skip a section
    read = 0

    # Make a new topology file for the system
    with open('topol_new.top','w') as newf:

        # Make a new .itp file for the monomer topology
    	with open(monomer+'.itp','w') as itpf:

            # Read the previous system topology line-by-line
    		with open('topol_mon.top','r') as oldf:
    			doc = oldf.readlines()
    			for line in doc:
    				if 'moleculetype' in line:
    					read = 1
    					newf.write('#include \"./'+monomer+'.itp\"\n\n')
    				if read == 0:
    					if 'Other' in line:
    						newf.write('Other         	'+n_monomers+'\n')
    					else:
    						newf.write(line)
    				elif '; Include Position restraint file' in line:
    					read = 0
    					newf.write(line)
    				else:
    					itpf.write(line)

    # Replace the topology file with the new one
    os.system('mv topol_new.top topol.top')
