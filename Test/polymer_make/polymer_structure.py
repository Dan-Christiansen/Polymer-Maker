import os, numpy as np
from pathlib import Path


def oss(fun):
    os.system(fun)


def polymer_structure(angle=120):
    monF = {}  # Keys are the characters, values are a tuple of the filenames and number of atoms
    with open('polymer_formula', 'r') as polyf:
        doc = polyf.readlines()
        # Polymer character string and DOP
        [formula, dop] = [doc[0].strip(), len(doc[0].strip())]
        for line in doc[1:]:
            tmp = line.split()
            this_dir, this_filename = os.path.split(__file__)
            data_path = Path(os.path.join(this_dir, "Molecule_Files/"))
            with open(str(data_path / (tmp[-1]+'.gro'))) as monf:
                mondoc = monf.readlines()
                atms = int(mondoc[1])
            monF[tmp[0]] = (tmp[1], atms)

    anumT = 0  # Total number of atoms in the system to be made
    # Input rotation angle between neighboring monomers
    rotr = 2*np.pi*float(angle)/360

    rotd = 0  # Start rotation is 0 degrees
    dx = 0.2
    xabs = 0
    natm = 0
    nmol = 0
    atmlist = []
    for mon in formula:
        anumT += monF[mon][1]
        nmol += 1
        xabs += dx
        with open(str(data_path / (monF[mon][0]+'.gro')), 'r') as grof:
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
            atmlist.append([str(nmol).rjust(
                5)+line[5:15]+str(natm).rjust(5), "%8.3f" % x, "%8.3f" % y, "%8.3f" % z])
        rotd += rotr
    xl = [float(i[1]) for i in atmlist]
    yl = [float(i[2]) for i in atmlist]
    zl = [float(i[3]) for i in atmlist]
    xlim = max(xl)-min(xl)
    ylim = max(yl)-min(yl)
    zlim = max(zl)-min(zl)
    with open('newmol.gro', 'w') as newf:
        newf.write('%d DOP\n%d\n' % (dop, anumT))
        for i in atmlist:
            newf.write('%s%s%s%s\n' % (i[0], i[1], i[2], i[3]))
        newf.write('%8.3f%8.3f%8.3f\n' % (xlim, ylim, zlim))
    oss('/usr/local/gromacs_gmx/bin/gmx editconf -f newmol.gro -box '+str(xlim)+' '+str(ylim)+' '+str(zlim)+' -c -o newmol.gro')
