# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:21:44 2021

@author: Dan
"""

import os

this_dir, this_filename = os.path.split(__file__)
data_path = os.path.join(this_dir, "Molecule_Files\\")
with open(data_path+"00240.gro",'r') as grof:
    doc = grof.readlines()
print(doc)