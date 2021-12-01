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

def polymer_clean():
    mol = 'newmol'

    with open(mol+'_new.itp','w') as newitpf:
    	with open(mol+'.itp','r') as itpf:
    		doc = itpf.readlines()
    	dihedralflag = 0
    	for line in doc:
    		if dihedralflag == 1:
    			if len(line.split()) >= 6:
    				newitpf.write(line)
    		else:
    			newitpf.write(line)
    		if '[ dihedrals ]' in line:
    			dihedralflag += 1
    os.system('mv '+mol+'_new.itp '+mol+'.itp')
