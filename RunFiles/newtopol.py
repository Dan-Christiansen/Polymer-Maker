import os, sys
mon = sys.argv[1] # Polymer filename, usually newmol
nmol = sys.argv[2] # Number of polymers that will be in the system
os.system('/usr/local/gromacs_gmx/bin/gmx pdb2gmx -f '+mon+'.gro -ff oplsaa -p topol_mon.top -water spce')
read = 0
with open('topol_new.top','w') as newf:
	with open(mon+'.itp','w') as itpf:
		with open('topol_mon.top','r') as oldf:
			doc = oldf.readlines()
			for line in doc:
				if 'moleculetype' in line:
					read = 1
					newf.write('#include \"./'+mon+'.itp\"\n\n')
				if read == 0:
					if 'Other' in line:
						newf.write('Other         	'+nmol+'\n')
					else:
						newf.write(line)
				elif '; Include Position restraint file' in line:
					read = 0
					newf.write(line)
				else:
					itpf.write(line)
os.system('mv topol_new.top topol.top')
#os.system('rm topol_mon.top')
