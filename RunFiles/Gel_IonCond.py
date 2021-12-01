import os, sys
from Full_Polymer_Make import polymer_clean, polymer_structure, polymer_topology

def oss(fun):
    os.system(fun)

def gmx(*fun):
	gmxrun = '/usr/local/gromacs_gmx/bin/gmx'
	for arg in fun:
		gmxrun += ' '+arg
	oss(gmxrun)

def main():
    # Prep the machine to use Gromacs
	oss('export PATH=$PATH:/usr/local/openmpi/bin')
	oss('export OMP_NUM_THREADS=4')

    with open('polymer_formula','w') as polymer_formula:
        polymer_formula.write('AAAAA\nA 00240')

    ## Need to figure out polymer structure generation
	with open('newmol.gro','r') as f: # Read polymer structure file for box size
		doc = f.readlines()
    box_x = doc[-1].split()[0] # Grab longest box length (always meant to be x-direction)
	gmx('editconf','-f newmol.gro','-box '+box_x+' 5 5','-c','-o box.gro') # Make the box the right size
	gmx('grompp','-f source/minim.mdp','-c box.gro','-p topol.top','-o em.tpr') # Minimize polymer structure to make backbone more correct
	gmx('mdrun','-deffnm em')
    gmx('insert-molecules','-ci em.gro','-nmol '+n_polymers,'-box 20 20 20','-o box.gro')
	gmx('solvate','-cp box.gro','-cs spc216.gro','-p topol.top','-o wetbox.gro') # Add SPCE water molecule (Chosen for its more accurate water dipole representation)
	gmx('grompp','-f source/minim.mdp','-c wetbox.gro','-p topol.top','-o prenvt.tpr')
	gmx('mdrun','-deffnm prenvt')
	gmx('grompp','-f source/npt.mdp','-c prenvt.gro','-p topol.top','-o npt.tpr')
	gmx('mdrun','-deffnm npt')
	nptflag = 0
	while not os.path.isfile('npt.gro'):
#		gmx('grompp','-f source/npt.mdp','-c em_ionwet.gro','-p topol.top','-o npt.tpr') # Perform a 1 ns NPT to equilibrate the system
		oss('rm \#*')
		gmx('grompp','-f source/npt.mdp','-c prenvt.gro','-p topol.top','-o npt.tpr') # Perform a 1 ns NPT to equilibrate the system
	        gmx('mdrun','-deffnm npt')
		nptflag += 1
		if nptflag == 10:
			break
	if conc != '0':
	######## CHANGED SECTION #########
		gmx('grompp','-f source/ion.mdp','-c npt.gro','-p topol.top','-o ions.tpr')
		oss('echo \"SOL\" | /usr/local/gromacs_gmx/bin/gmx genion -s ions.tpr -pname '+catType+' -nname '+anType+' -conc '+conc+' -o ionwetbox.gro -p topol.top') # Add desired ions

	nvtflag = 0
	while not os.path.isfile('nvt.gro'):
		if conc != '0':
			gmx('grompp','-f source/nvt.mdp','-c ionwetbox.gro','-p topol.top','-o nvt.tpr') # Production run 10 ns NVT
		else:
			gmx('grompp','-f source/nvt.mdp','-c npt.gro','-p topol.top','-o nvt.tpr')
		gmx('mdrun','-deffnm nvt')
		nvtflag += 1
		if nvtflag == 10:
			break

if __name__ == "__main__":
    main()
