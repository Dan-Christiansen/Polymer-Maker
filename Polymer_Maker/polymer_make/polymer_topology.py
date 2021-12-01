import os, numpy as np
from pathlib import Path


def oss(fun):  # Abbreviated 'os.system' function
    os.system(fun)


def polymer_topology(n_monomers):
    monomer = 'newmol'  # Polymer filename, usually newmol

    this_dir, this_filename = os.path.split(__file__)
    os.environ['GMXLIB'] = str(this_dir)
#    os.system('/usr/local/gromacs_gmx/bin/gmx pdb2gmx -f '+monomer+'.gro -ff '+str(data_path)+' -p topol_mon.top -water spce')
    os.system('echo \"1\" | /usr/local/gromacs_gmx/bin/gmx pdb2gmx -f '+monomer+'.gro -p topol_mon.top -water spce')
#    oss('unset GMXLIB')
    del os.environ['GMXLIB']

    # Flag for when to skip a section
    read = 0

    # Make a new topology file for the system
    with open('topol_new.top', 'w') as newf:

        # Make a new .itp file for the monomer topology
        with open(monomer+'.itp', 'w') as itpf:

            # Read the previous system topology line-by-line
            with open('topol_mon.top', 'r') as oldf:
                doc = oldf.readlines()
                for line in doc:
                    if 'moleculetype' in line:
                        read = 1
                        newf.write('#include \"./'+monomer+'.itp\"\n\n')
                    if read == 0:
                        if 'Other' in line:
                            newf.write('Other         	'+str(n_monomers)+'\n')
                        else:
                            newf.write(line)
                    elif '; Include Position restraint file' in line:
                        read = 0
                        newf.write(line)
                    else:
                        itpf.write(line)

    # Replace the topology file with the new one
    os.system('mv topol_new.top topol.top')
