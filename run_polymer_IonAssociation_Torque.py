from collections import Counter as cnt
import os
cations = ['LI','K']
anions = ['CL']
concentrations = ['0','0.2','0.05','0.1']
r = 5
molnames = list(set([i.split('.')[0] for i in os.listdir(os.getcwd()+'/Molecule_Files/') if '.gro' in i]))

job = 1
#for f in molnames:
#	for conc in concentrations:
#		for cat in cations:
#			for an in anions:
for cat in cations:
	for an in anions:
		for conc in concentrations:
			for f in molnames:
				for trial in range(1,r+1):
					outname = f+'_'+cat+'_'+an+'_'+conc.replace('.','p')+'M_'+str(trial)
					os.mkdir(outname)
					os.system('cp -r RunFiles/* '+outname)
					os.system('cp -r ./Molecule_Files/'+f+'.* '+outname+'/source/')
					os.chdir(outname)
					with open('polymer_formula','w') as pff:
						pff.write('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nA '+f+'\n')

					with open('dag.dag','w') as dagf:
						job += 1
						dagf.write('  JOB  job'+str(job)+'  submit.sub\n')
						dagf.write('  VARS job'+str(job)+' molecule=\"'+f+'\" cation=\"'+cat+'\" anion=\"'+an+'\" conc=\"'+conc+'\" r=\"'+str(trial)+'\" outname=\"'+outname+'\"\n')
					os.system('condor_submit_dag dag.dag')
					os.chdir('..')
#					exit()
