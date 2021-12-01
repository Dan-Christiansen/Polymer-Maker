#!/usr/bin/python
#############################################################################
### Ion association of zwitterionic polymers
### Developed by Daniel Christiansen
### Feb. 17th, 2021
###
### Intended to be run on HTCONDOR
#############################################################################
import os, sys, tarfile, shutil
from collections import Counter as cnt

# Dictionary of monomer types and atoms of interest (i.e., atoms that will be used as reference for the radial distribution function, coordination number, etc.)
# Each key has an array as a key. This array consists of a value and two subarrays. The value is the total number of atoms in the monomer.
# The first array belongs to the positively charged zw. moiety, the latter to the negatively charged.
# This will absolutely need to be updated when new monomers are used.
#zwdict = {'CBAA':[33,[20,24],[32,33]],'SBMA':[39,[19,23],[37,38,39]]}
#zwdict = {'CBAA1':[33,[20,24],[32,33]],'CBAA2':[36,[20,24],[35,36]],'SBAA1':[34,[20,24],[32,33,34]],'SBAA2':[37,[20,24],[35,36,37]]}
zwdict = {'00412': [['C09'], ['O0F', 'O0G', 'O0H']], '00413': [['C09'], ['O0H', 'O0I', 'O0J']], '00410': [['C09'], ['O0M', 'O0N', 'O0O']], '00411': [['C09'], ['O0D', 'O0E', 'O0F']], '00417': [['C09'], ['O0D', 'O0E', 'O0F']], '00414': [['C09'], ['O0J', 'O0K', 'O0M']], '00415': [['C09'], ['O0M', 'O0N', 'O0O']], '00418': [['C09'], ['O0E', 'O0F', 'O0G']], '00419': [['C09'], ['O0F', 'O0G', 'O0H']], '00502': [['C0G', 'C0H'], ['O0A', 'O0B']], '00503': [['C0I', 'C0J'], ['O0A', 'O0B']], '00319': [['H09', 'C0A'], ['O0N', 'O0O', 'O0P']], '00318': [['H09', 'C0A'], ['O0K', 'O0M', 'O0N']], '00314': [['H09', 'C0A'], ['O0N', 'O0O', 'O0P']], '00317': [['H09', 'C0A'], ['O0I', 'O0J', 'O0K']], '00316': [['H09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00311': [['H09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00310': [['H09', 'C0A'], ['O0E', 'O0F', 'O0G']], '00313': [['H09', 'C0A'], ['O0K', 'O0M', 'O0N']], '00312': [['H09', 'C0A'], ['O0I', 'O0J', 'O0K']], '00395': [['C09'], ['O0D', 'O0E']], '00394': [['C09'], ['O0M', 'O0N']], '00397': [['C09'], ['O0H', 'O0I']], '00396': [['C09'], ['O0F', 'O0G']], '00258': [['C09', 'C0A'], ['O0E', 'O0F', 'O0G']], '00259': [['C09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00393': [['C09'], ['O0J', 'O0K']], '00392': [['C09'], ['O0H', 'O0I']], '00254': [['C09', 'C0A'], ['O0K', 'O0M']], '00255': [['C09', 'C0A'], ['O0N', 'O0O']], '00490': [['H0P', 'H0Q', 'H0R'], ['O0A', 'O0B']], '00257': [['C09', 'C0A'], ['O0D', 'O0E', 'O0F']], '00250': [['C09', 'C0A'], ['O0N', 'O0O']], '00251': [['C09', 'C0A'], ['O0E', 'O0F']], '00252': [['C09', 'C0A'], ['O0G', 'O0H']], '00253': [['C09', 'C0A'], ['O0I', 'O0J']], '00429': [['C09'], ['O0I', 'O0J', 'O0K']], '00427': [['C09'], ['O0E', 'O0F', 'O0G']], '00426': [['C09'], ['O0N', 'O0O', 'O0P']], '00425': [['C09'], ['O0K', 'O0M', 'O0N']], '00424': [['C09'], ['O0I', 'O0J', 'O0K']], '00423': [['C09'], ['O0G', 'O0H', 'O0I']], '00422': [['C09'], ['O0E', 'O0F', 'O0G']], '00421': [['C09'], ['O0H', 'O0I', 'O0J']], '00420': [['C09'], ['O0G', 'O0H', 'O0I']], '00351': [['H09', 'H0A'], ['O0N', 'O0O']], '00308': [['H09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00309': [['H09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00306': [['H09', 'C0A'], ['O0E', 'O0F', 'O0G']], '00307': [['H09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00304': [['H09', 'C0A'], ['O0C', 'O0D', 'O0E']], '00305': [['H09', 'C0A'], ['O0D', 'O0E', 'O0F']], '00302': [['H09', 'C0A'], ['O0K', 'O0M']], '00303': [['H09', 'C0A'], ['O0N', 'O0O']], '00300': [['H09', 'C0A'], ['O0G', 'O0H']], '00301': [['H09', 'C0A'], ['O0I', 'O0J']], '00386': [['C09'], ['O0D', 'O0E']], '00387': [['C09'], ['O0E', 'O0F']], '00384': [['C09'], ['O0B', 'O0C']], '00385': [['C09'], ['O0C', 'O0D']], '00382': [['H09', 'H0A'], ['O0M', 'O0N', 'O0O']], '00383': [['H09', 'H0A'], ['O0O', 'O0P', 'O0Q']], '00380': [['H09', 'H0A'], ['O0H', 'O0I', 'O0J']], '00381': [['H09', 'H0A'], ['O0J', 'O0K', 'O0M']], '00350': [['H09', 'H0A'], ['O0K', 'O0M']], '00388': [['C09'], ['O0F', 'O0G']], '00389': [['C09'], ['O0G', 'O0H']], '00451': [['H0H', 'C0I', 'C0J'], ['O0A', 'O0B']], '00496': [['C0E', 'C0F'], ['O0A', 'O0B']], '00409': [['C09'], ['O0J', 'O0K', 'O0M']], '00359': [['H09', 'H0A'], ['O0G', 'O0H', 'O0I']], '00498': [['C0G', 'C0H'], ['O0A', 'O0B']], '00358': [['H09', 'H0A'], ['O0E', 'O0F', 'O0G']], '00430': [['C09'], ['O0K', 'O0M', 'O0N']], '00431': [['C09'], ['O0N', 'O0O', 'O0P']], '00432': [['C0E', 'C0F', 'C0G'], ['O0A', 'O0B']], '00433': [['C0F', 'C0G', 'C0H'], ['O0A', 'O0B']], '00434': [['C0G', 'C0H', 'C0I'], ['O0A', 'O0B']], '00435': [['C0H', 'C0I', 'C0J'], ['O0A', 'O0B']], '00436': [['C0I', 'C0J', 'C0K'], ['O0A', 'O0B']], '00437': [['C0J', 'C0K', 'C0M'], ['O0A', 'O0B']], '00339': [['H09', 'H0A'], ['O0F', 'O0G']], '00338': [['H09', 'H0A'], ['O0E', 'O0F']], '00333': [['H09', 'C0A'], ['O0J', 'O0K', 'O0M']], '00332': [['H09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00331': [['H09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00330': [['H09', 'C0A'], ['O0O', 'O0P', 'O0Q']], '00337': [['H09', 'H0A'], ['O0D', 'O0E']], '00336': [['H09', 'H0A'], ['O0C', 'O0D']], '00335': [['H09', 'C0A'], ['O0O', 'O0P', 'O0Q']], '00334': [['H09', 'C0A'], ['O0M', 'O0N', 'O0O']], '00511': [['C0P', 'C0Q'], ['O0A', 'O0B']], '00510': [['C0N', 'C0O'], ['O0A', 'O0B']], '00324': [['H09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00325': [['H09', 'C0A'], ['O0I', 'O0J', 'O0K']], '00326': [['H09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00327': [['H09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00320': [['H09', 'C0A'], ['O0D', 'O0E', 'O0F']], '00321': [['H09', 'C0A'], ['O0E', 'O0F', 'O0G']], '00322': [['H09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00288': [['H09', 'C0A'], ['O0C', 'O0D']], '00287': [['C09', 'C0A'], ['O0O', 'O0P', 'O0Q']], '00286': [['C09', 'C0A'], ['O0M', 'O0N', 'O0O']], '00285': [['C09', 'C0A'], ['O0J', 'O0K', 'O0M']], '00284': [['C09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00328': [['H09', 'C0A'], ['O0J', 'O0K', 'O0M']], '00329': [['H09', 'C0A'], ['O0M', 'O0N', 'O0O']], '00281': [['C09', 'C0A'], ['O0M', 'O0N', 'O0O']], '00280': [['C09', 'C0A'], ['O0J', 'O0K', 'O0M']], '00400': [['C09'], ['O0B', 'O0C', 'O0D']], '00445': [['C0K', 'C0M', 'C0N'], ['O0A', 'O0B']], '00444': [['C0I', 'C0J', 'C0K'], ['O0A', 'O0B']], '00497': [['C0F', 'C0G'], ['O0A', 'O0B']], '00446': [['C0N', 'C0O', 'C0P'], ['O0A', 'O0B']], '00441': [['C0N', 'C0O', 'C0P'], ['O0A', 'O0B']], '00440': [['C0K', 'C0M', 'C0N'], ['O0A', 'O0B']], '00457': [['H0N', 'C0O', 'C0P'], ['O0A', 'O0B']], '00289': [['H09', 'C0A'], ['O0D', 'O0E']], '00323': [['H09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00298': [['H09', 'C0A'], ['O0N', 'O0O']], '00299': [['H09', 'C0A'], ['O0E', 'O0F']], '00353': [['H09', 'H0A'], ['O0D', 'O0E', 'O0F']], '00352': [['H09', 'H0A'], ['O0C', 'O0D', 'O0E']], '00355': [['H09', 'H0A'], ['O0F', 'O0G', 'O0H']], '00354': [['H09', 'H0A'], ['O0E', 'O0F', 'O0G']], '00357': [['H09', 'H0A'], ['O0H', 'O0I', 'O0J']], '00356': [['H09', 'H0A'], ['O0G', 'O0H', 'O0I']], '00290': [['H09', 'C0A'], ['O0E', 'O0F']], '00291': [['H09', 'C0A'], ['O0F', 'O0G']], '00292': [['H09', 'C0A'], ['O0G', 'O0H']], '00293': [['H09', 'C0A'], ['O0H', 'O0I']], '00294': [['H09', 'C0A'], ['O0E', 'O0F']], '00295': [['H09', 'C0A'], ['O0G', 'O0H']], '00296': [['H09', 'C0A'], ['O0I', 'O0J']], '00297': [['H09', 'C0A'], ['O0K', 'O0M']], '00458': [['H0P', 'C0Q', 'C0R'], ['O0A', 'O0B']], '00459': [['H0G', 'C0H', 'C0I'], ['O0A', 'O0B']], '00456': [['H0K', 'C0M', 'C0N'], ['O0A', 'O0B']], '00283': [['C09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00455': [['H0I', 'C0J', 'C0K'], ['O0A', 'O0B']], '00452': [['H0I', 'C0J', 'C0K'], ['O0A', 'O0B']], '00453': [['H0J', 'C0K', 'C0M'], ['O0A', 'O0B']], '00450': [['H0G', 'C0H', 'C0I'], ['O0A', 'O0B']], '00282': [['C09', 'C0A'], ['O0O', 'O0P', 'O0Q']], '00499': [['C0H', 'C0I'], ['O0A', 'O0B']], '00391': [['C09'], ['O0F', 'O0G']], '00390': [['C09'], ['O0D', 'O0E']], '00507': [['C0G', 'C0H'], ['O0A', 'O0B']], '00492': [['H0I', 'H0J', 'H0K'], ['O0A', 'O0B']], '00342': [['H09', 'H0A'], ['O0E', 'O0F']], '00343': [['H09', 'H0A'], ['O0G', 'O0H']], '00340': [['H09', 'H0A'], ['O0G', 'O0H']], '00341': [['H09', 'H0A'], ['O0H', 'O0I']], '00346': [['H09', 'H0A'], ['O0N', 'O0O']], '00347': [['H09', 'H0A'], ['O0E', 'O0F']], '00344': [['H09', 'H0A'], ['O0I', 'O0J']], '00345': [['H09', 'H0A'], ['O0K', 'O0M']], '00500': [['C0I', 'C0J'], ['O0A', 'O0B']], '00501': [['C0J', 'C0K'], ['O0A', 'O0B']], '00348': [['H09', 'H0A'], ['O0G', 'O0H']], '00349': [['H09', 'H0A'], ['O0I', 'O0J']], '00504': [['C0K', 'C0M'], ['O0A', 'O0B']], '00505': [['C0N', 'C0O'], ['O0A', 'O0B']], '00506': [['C0P', 'C0Q'], ['O0A', 'O0B']], '00462': [['H0N', 'C0O', 'C0P'], ['O0A', 'O0B']], '00469': [['H0J', 'H0K', 'C0M'], ['O0A', 'O0B']], '00468': [['H0I', 'H0J', 'C0K'], ['O0A', 'O0B']], '00493': [['H0K', 'H0M', 'H0N'], ['O0A', 'O0B']], '00439': [['C0I', 'C0J', 'C0K'], ['O0A', 'O0B']], '00463': [['H0P', 'C0Q', 'C0R'], ['O0A', 'O0B']], '00399': [['C09'], ['O0M', 'O0N']], '00467': [['H0H', 'H0I', 'C0J'], ['O0A', 'O0B']], '00466': [['H0G', 'H0H', 'C0I'], ['O0A', 'O0B']], '00465': [['H0F', 'H0G', 'C0H'], ['O0A', 'O0B']], '00398': [['C09'], ['O0J', 'O0K']], '00265': [['C09', 'C0A'], ['O0K', 'O0M', 'O0N']], '00264': [['C09', 'C0A'], ['O0I', 'O0J', 'O0K']], '00267': [['C09', 'C0A'], ['O0E', 'O0F', 'O0G']], '00266': [['C09', 'C0A'], ['O0N', 'O0O', 'O0P']], '00261': [['C09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00260': [['C09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00263': [['C09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00262': [['C09', 'C0A'], ['O0E', 'O0F', 'O0G']], '00269': [['C09', 'C0A'], ['O0I', 'O0J', 'O0K']], '00268': [['C09', 'C0A'], ['O0G', 'O0H', 'O0I']], '00408': [['C09'], ['O0H', 'O0I', 'O0J']], '00472': [['H0K', 'H0M', 'C0N'], ['O0A', 'O0B']], '00483': [['H0H', 'H0I', 'H0J'], ['O0A', 'O0B']], '00473': [['H0N', 'H0O', 'C0P'], ['O0A', 'O0B']], '00494': [['H0N', 'H0O', 'H0P'], ['O0A', 'O0B']], '00474': [['H0P', 'H0Q', 'C0R'], ['O0A', 'O0B']], '00475': [['H0G', 'H0H', 'C0I'], ['O0A', 'O0B']], '00476': [['H0I', 'H0J', 'C0K'], ['O0A', 'O0B']], '00477': [['H0K', 'H0M', 'C0N'], ['O0A', 'O0B']], '00470': [['H0G', 'H0H', 'C0I'], ['O0A', 'O0B']], '00471': [['H0I', 'H0J', 'C0K'], ['O0A', 'O0B']], '00379': [['H09', 'H0A'], ['O0F', 'O0G', 'O0H']], '00378': [['H09', 'H0A'], ['O0O', 'O0P', 'O0Q']], '00377': [['H09', 'H0A'], ['O0M', 'O0N', 'O0O']], '00376': [['H09', 'H0A'], ['O0J', 'O0K', 'O0M']], '00375': [['H09', 'H0A'], ['O0H', 'O0I', 'O0J']], '00374': [['H09', 'H0A'], ['O0F', 'O0G', 'O0H']], '00373': [['H09', 'H0A'], ['O0I', 'O0J', 'O0K']], '00372': [['H09', 'H0A'], ['O0H', 'O0I', 'O0J']], '00371': [['H09', 'H0A'], ['O0G', 'O0H', 'O0I']], '00370': [['H09', 'H0A'], ['O0F', 'O0G', 'O0H']], '00489': [['H0N', 'H0O', 'H0P'], ['O0A', 'O0B']], '00488': [['H0K', 'H0M', 'H0N'], ['O0A', 'O0B']], '00478': [['H0N', 'H0O', 'C0P'], ['O0A', 'O0B']], '00479': [['H0P', 'H0Q', 'C0R'], ['O0A', 'O0B']], '00276': [['C09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00277': [['C09', 'C0A'], ['O0I', 'O0J', 'O0K']], '00270': [['C09', 'C0A'], ['O0K', 'O0M', 'O0N']], '00278': [['C09', 'C0A'], ['O0F', 'O0G', 'O0H']], '00279': [['C09', 'C0A'], ['O0H', 'O0I', 'O0J']], '00491': [['H0G', 'H0H', 'H0I'], ['O0A', 'O0B']], '00401': [['C09'], ['O0C', 'O0D', 'O0E']], '00369': [['H09', 'H0A'], ['O0E', 'O0F', 'O0G']], '00402': [['C09'], ['O0D', 'O0E', 'O0F']], '00405': [['C09'], ['O0G', 'O0H', 'O0I']], '00404': [['C09'], ['O0F', 'O0G', 'O0H']], '00407': [['C09'], ['O0F', 'O0G', 'O0H']], '00406': [['C09'], ['O0D', 'O0E', 'O0F']], '00360': [['H09', 'H0A'], ['O0I', 'O0J', 'O0K']], '00361': [['H09', 'H0A'], ['O0K', 'O0M', 'O0N']], '00362': [['H09', 'H0A'], ['O0N', 'O0O', 'O0P']], '00363': [['H09', 'H0A'], ['O0E', 'O0F', 'O0G']], '00364': [['H09', 'H0A'], ['O0G', 'O0H', 'O0I']], '00365': [['H09', 'H0A'], ['O0I', 'O0J', 'O0K']], '00366': [['H09', 'H0A'], ['O0K', 'O0M', 'O0N']], '00367': [['H09', 'H0A'], ['O0N', 'O0O', 'O0P']], '00508': [['C0I', 'C0J'], ['O0A', 'O0B']], '00482': [['H0G', 'H0H', 'H0I'], ['O0A', 'O0B']], '00509': [['C0K', 'C0M'], ['O0A', 'O0B']], '00243': [['C09', 'C0A'], ['O0F', 'O0G']], '00242': [['C09', 'C0A'], ['O0E', 'O0F']], '00241': [['C09', 'C0A'], ['O0D', 'O0E']], '00240': [['C09', 'C0A'], ['O0C', 'O0D']], '00247': [['C09', 'C0A'], ['O0G', 'O0H']], '00246': [['C09', 'C0A'], ['O0E', 'O0F']], '00245': [['C09', 'C0A'], ['O0H', 'O0I']], '00244': [['C09', 'C0A'], ['O0G', 'O0H']], '00481': [['H0F', 'H0G', 'H0H'], ['O0A', 'O0B']], '00249': [['C09', 'C0A'], ['O0K', 'O0M']], '00248': [['C09', 'C0A'], ['O0I', 'O0J']], '00485': [['H0J', 'H0K', 'H0M'], ['O0A', 'O0B']], '00484': [['H0I', 'H0J', 'H0K'], ['O0A', 'O0B']], '00487': [['H0I', 'H0J', 'H0K'], ['O0A', 'O0B']], '00486': [['H0G', 'H0H', 'H0I'], ['O0A', 'O0B']],'00256':[['C09','C0A'],['O0C','O0D','O0E']],'00271':[['C09','C0A'],['O0N','O0O','O0P']],'00272':[['C09','C0A'],['O0D','O0E','O0F']],'00273':[['C09','C0A'],['O0E','O0F','O0G']],'00274':[['C09','C0A'],['O0F','O0G','O0H']],'00275':[['C09','C0A'],['O0G','O0H','O0I']],'00315':[['H09','C0A'],['O0E','O0F','O0G']],'00368':[['H09','H0A'],['O0D','O0E','O0F']],'00403':[['C09'],['O0E','O0F','O0G']],'00416':[['C09'],['O0C','O0D','O0E']],'00428':[['C09'],['O0G','O0H','O0I']],'00438':[['O0A','O0B'],['C0G','C0H','C0I']],'00442':[['O0A','O0B'],['C0P','C0Q','C0R']],'00443':[['O0A','O0B'],['C0G','C0H','C0I']],'00447':[['O0A','O0B'],['C0P','C0Q','C0R']],'00448':[['O0A','O0B'],['H0E','C0F','C0G']],'00449':[['O0A','O0B'],['H0F','C0G','C0H']],'00454':[['O0A','O0B'],['H0G','C0H','C0I']],'00460':[['O0A','O0B'],['H0I','C0J','C0K']],'00461':[['O0A','O0B'],['H0K','C0M','C0N']],'00464':[['O0A','O0B'],['H0E','H0F','C0G']],'00480':[['O0A','O0B'],['H0E','H0F','H0G']],'00495':[['O0A','O0B'],['H0P','H0Q','H0R']]}

# Execute bash commands with os.system function
def oss(fun):
	os.system(fun)

# Shortcut to make Gromacs commands more similar to style used in terminal (In my opinion)
# Ideally this would be replaced by the gmx python module and used as that is meant to be
def gmx(*fun):
	gmxrun = '/usr/local/gromacs_gmx/bin/gmx'
	for arg in fun:
		gmxrun += ' '+arg
	oss(gmxrun)

if __name__ == "__main__":
	# Prep the machine to use Gromacs
	oss('export PATH=$PATH:/usr/local/openmpi/bin')
	oss('export OMP_NUM_THREADS=4')

	# Ion types must match ions.itp 'molname' format in oplsaa.ff folder
	catType = sys.argv[1] # Cation type submitted as a string
	anType = sys.argv[2] # Anion type submitted as a string
	conc = sys.argv[3] # A float value submitted as a string defining the concentration of ions
	r = sys.argv[4] # Trial number
	oss('rm \#*')
	oss('python3 polymake_v3.py 120') # Create polymer structure
	oss('python newtopol.py newmol 1') # Create polymer .itp file and topol.top topology file (Note: OPLS-AA forcefield and SPCE water model chosen)
	oss('python remove_missing_dihedrals.py newmol') # Remove illegal dihedrals produced during polymerization. These dihedrals are the backbone carbons across any 3 neighboring monomers
	with open('newmol.gro','r') as f: # Read polymer structure file for box size
		doc = f.readlines()
	box_x = doc[-1].split()[0] # Grab longest box length (always meant to be x-direction)
	gmx('editconf','-f newmol.gro','-box '+box_x+' 5 5','-c','-o box.gro') # Make the box the right size
	gmx('grompp','-f source/minim.mdp','-c box.gro','-p topol.top','-o em.tpr') # Minimize polymer structure to make backbone more correct
	gmx('mdrun','-deffnm em')
	gmx('solvate','-cp em.gro','-cs spc216.gro','-p topol.top','-o wetbox.gro') # Add SPCE water molecule (Chosen for its more accurate water dipole representation)
## Adding ion before NPT - This may have caused some issues with ions beginning the NVT already associated
#	gmx('grompp','-f source/ion.mdp','-c wetbox.gro','-p topol.top','-o ions.tpr')
#	oss('echo \"SOL\" | /usr/local/gromacs_gmx/bin/gmx genion -s ions.tpr -pname '+catType+' -nname '+anType+' -conc '+conc+' -o ionwetbox.gro -p topol.top') # Add desired ions
#	gmx('grompp','-f source/minim.mdp','-c ionwetbox.gro','-p topol.top','-o em_ionwet.tpr') # Minimize system energy again now that it has water and ions
#	gmx('mdrun','-deffnm em_ionwet')
#	gmx('grompp','-f source/npt.mdp','-c em_ionwet.gro','-p topol.top','-o npt.tpr') # Perform a 1 ns NPT to equilibrate the system
######## CHANGED SECTION #########
	gmx('grompp','-f source/minim.mdp','-c wetbox.gro','-p topol.top','-o prenvt.tpr')
	gmx('mdrun','-deffnm prenvt')
	gmx('grompp','-f source/npt.mdp','-c prenvt.gro','-p topol.top','-o npt.tpr')
######## CHANGED SECTION #########
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

	with open('polymer_formula','r') as f:
		doc = f.readlines()
	formula = doc[0].strip('\n') # Polymer formula encoded
	monkeys = {i.split()[0]:i.split()[1] for i in doc[1:]} # Cypher for polymer formula relating letter codes to monomer filename

	index = 'printf \"'
	zwlist = []
	cnt = 0
	for i in formula:
		cnt += 1
		atoms = zwdict[monkeys[i]]
		cats = 'ri '+str(cnt)+' & a '
		catlist = 'r_'+str(cnt)
		for cat in atoms[0]:
			cats += cat+' '
			catlist += '_'+cat
		zwlist.append(catlist)
		index += cats + ' \n '
		ans = 'ri '+str(cnt)+' & a '
		anlist = 'r_'+str(cnt)
		for an in atoms[1]:
			ans += an+' '
			anlist += '_'+an
		zwlist.append(anlist)
		index += ans + ' \n '
	index += 'r '+catType+' to_use \n r '+anType+' to_use \n q \n\"'
	oss(index+' | /usr/local/gromacs_gmx/bin/gmx make_ndx -f nvt.gro')
	with open('index.ndx','r') as idxf:
		doc = idxf.readlines()
	with open('index.ndx','w') as idxf:
		for line in doc:
			if '&' not in line:
				idxf.write(line)
			else:
				idxf.write(line.split('&_')[0]+line.split('&_')[1])
	for z in zwlist:
		gmx('rdf','-f nvt.trr','-s nvt.gro','-n index.ndx','-o rdf_'+z+'_Ow.xvg','-cn rdf_'+z+'_Ow_CN.xvg','-ref '+z,'-sel SOL')
		if conc != '0':
			gmx('rdf','-f nvt.trr','-s nvt.gro','-n index.ndx','-o rdf_'+z+'_'+anType+'.xvg','-cn '+z+'_'+anType+'_CN.xvg','-ref '+z,'-sel '+anType+'_TO_USE')
			gmx('rdf','-f nvt.trr','-s nvt.gro','-n index.ndx','-o rdf_'+z+'_'+catType+'.xvg','-cn '+z+'_'+catType+'_CN.xvg','-ref '+z,'-sel '+catType+'_TO_USE')
	if conc != '0':
		os.system('echo \"'+anType+'_TO_USE\" | gmx msd -f nvt.trr -s nvt.gro -n index.ndx -o msd_'+anType+'.xvg')
		os.system('echo \"'+catType+'_TO_USE\" | gmx msd -f nvt.trr -s nvt.gro -n index.ndx -o msd_'+catType+'.xvg')

	# Files to keep
	keep = ['polymer_formula','npt.edr','nvt.gro','nvt.trr','rdf_zwcation_Ow.xvg','rdf_zwcation_'+anType+'.xvg','rdf_zwanion_Ow.xvg','rdf_zwanion_'+catType+'.xvg','index.ndx']

	files = os.listdir(os.getcwd())
	for f in files:
		if f in keep or '.xvg' in f:
			pass
		else:
			if os.path.isdir(f):
				shutil.rmtree(f,ignore_errors=True)
			else:
				os.remove(f)

	# Make a compressed .tgz file with desired files
#	tfname = 'ionassoc'
#	for i in range(len(cnt(formula).keys())):
#		tfname += '_'+cnt(formula).keys()[i]+'_'+str(cnt(formula).values()[i])
##	tf = tarfile.open(tfname+'_'+catType+'_'+anType+'_'+conc+'_'+r+'.tgz', "w:gz")
#	oss('mkdir logouterr')
#	oss('mv *err *out logouterr')
#	tf.add('logouterr')
#	oss('rm -r logouterr oplsaa.ff source')
#	files = os.listdir(os.getcwd())
#	files.remove(tfname+'_'+catType+'_'+anType+'_'+conc+'_'+r+'.tgz')
#	for f in files:
#		tf.add(f)
#		if f in keep:
#			tf.add(f)
#			os.remove(f)
#		else:
#			os.remove(f)
#	tf.close()
