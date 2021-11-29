########################################################################################
#
# 				  Polymer Maker v.3
#
#		Make a linear polymer with a more accurate rotation of 
#		the monomers about the backbone. This is still not a
#		'realistic' method, where ab initio would be more 
#		approriate. However, it's better than having all monomers
#		aligned in one direction.
#
#			- DC
#
########################################################################################
# 			     Update from v.2 (Feb. 4th 2021)
# - Include definable monomer type and order from an input file called 'copolymer_formula'
# - Controllable rotation parameter
#
# Format of the input file, 'polymer_formula':
#	- First line: String of single characters, where each character corresponds to
#		      a single monomer
#	- Other lines: Every line after the first consists of two columns. The first 
#		       column is the single character used to identify the monomer;
#		       the second is the corresponding monomer filename without '.gro'
########################################################################################

import os, sys
import numpy as np

def oss(fun): # Abbreviated 'os.system' function
	os.system(fun)

def mainfunc():
	oss('rm \#*') # Remove temporary files if there are any

	# Read the formula and filenames of the system from 'polymer_formula' file.
	# Creates a dictionary of the monomer characters, filenames, and atoms per monomer
	monF = {} # Keys are the characters, values are a tuple of the filenames and number of atoms
	with open('polymer_formula','r') as polyf:
		doc = polyf.readlines()
		[formula, dop] = [doc[0].strip(),len(doc[0].strip())] # Polymer character string and DOP
		for line in doc[1:]:
			tmp = line.split()
			with open('./source/'+tmp[-1]+'.gro') as monf:
				mondoc = monf.readlines()
				atms = int(mondoc[1])
			monF[tmp[0]] = (tmp[1],atms)

	anumT = 0 # Total number of atoms in the system to be made
	rotr = float(sys.argv[1]) # Input rotation
	rotr = 2*np.pi*rotr/360

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
		with open('./source/'+monF[mon][0]+'.gro','r') as grof:
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
if __name__ == '__main__':
	mainfunc()
