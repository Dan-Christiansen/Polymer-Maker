# Tool to make linear polymer structure and topology files for use with GROMACS OPLS-AA molecular dynamics simulation
#
#
#
# Input files: polymer_formula
# Input arguments:
# Number of polymers in the system
# Rotation angle between neighboring monomers
#
import os, sys
import numpy as np

def oss(fun): # Abbreviated 'os.system' function
    os.system(fun)

def polymer_structure():
    # Read the formula and filenames of the system from 'polymer_formula' file.
    # Creates a dictionary of the monomer characters, filenames, and atoms per monomer
    monF = {} # Keys are the characters, values are a tuple of the filenames and number of atoms
    with open('polymer_formula','r') as polyf:
        doc = polyf.readlines()
        [formula, dop] = [doc[0].strip(),len(doc[0].strip())] # Polymer character string and DOP
        for line in doc[1:]:
            tmp = line.split()
            with open('./Molecule_Files/'+tmp[-1]+'.gro') as monf:
                mondoc = monf.readlines()
                atms = int(mondoc[1])
            monF[tmp[0]] = (tmp[1],atms)

    anumT = 0 # Total number of atoms in the system to be made
    rotr = 2*np.pi*float(sys.argv[2])/360 # Input rotation angle between neighboring monomers

    rotd = 0 # Start rotation is 0 degrees
    dx = 0.2
    xabs = 0
    natm = 0
    nmol = 0
    atmlist = []
    for mon in formula:
        anumT += monF[mon][1]
        nmol += 1
        xabs += dx
        with open('./Molecule_Files/'+monF[mon][0]+'.gro','r') as grof:
            doc = grof.readlines()
        for line in doc[2:-1]:
            if 'C01' in line:
                pos = [float(i) for i in line.split()[3:]]
                yshift = -1*(pos[1]*np.cos(rotr)-pos[2]*np.sin(rotr))
                zshift = -1*(pos[1]*np.sin(rotr)+pos[2]*np.cos(rotr))
                break
        for line in doc[2:-1]:
            natm += 1
            pos = [float(i) for i in line.split()[3:]]
            x = pos[0]+xabs
            y = pos[1]*np.cos(rotd)-pos[2]*np.sin(rotd)+yshift
            z = pos[1]*np.sin(rotd)+pos[2]*np.cos(rotd)+zshift
            #newf.write('%s%8.3f%8.3f%8.3f\n' % (line[:20],pos[0],y,z))
            atmlist.append([str(nmol).rjust(5)+line[5:15]+str(natm).rjust(5),"%8.3f" % x,"%8.3f" % y,"%8.3f" % z])
        rotd += rotr
    xl = [float(i[1]) for i in atmlist]
    yl = [float(i[2]) for i in atmlist]
    zl = [float(i[3]) for i in atmlist]
    xlim = max(xl)-min(xl)
    ylim = max(yl)-min(yl)
    zlim = max(zl)-min(zl)
    with open('newmol.gro','w') as newf:
        newf.write('%d DOP\n%d\n' % (dop,anumT))
        for i in atmlist:
            newf.write('%s%s%s%s\n' % (i[0],i[1],i[2],i[3]))
        newf.write('%8.3f%8.3f%8.3f\n' % (xlim,ylim,zlim))
    oss('/usr/local/gromacs_gmx/bin/gmx editconf -f newmol.gro -box '+str(xlim)+' '+str(ylim)+' '+str(zlim)+' -c -o newmol.gro')

def polymer_topology():
    monomer = 'newmol' # Polymer filename, usually newmol
    n_monomers = sys.argv[1] # Number of polymers that will be in the system

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


def polymer_make():
    oss('rm \#*') # Remove temporary files if there are any
    polymer_structure()
    polymer_topology()
    polymer_clean()

if __name__ == '__main__':
    polymer_make()
