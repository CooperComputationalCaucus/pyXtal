'''
Script to pull in set of csv files and make on single dataframe .

@author: pmm
'''
import numpy as np
import os
import sys
sys.path.append('../')
import pyXtal as pxt
from pyXtal.csp_utils.dataset_io import csv_to_dataframe

func = pxt.csp_utils.dataset_io.replace_spaces
dir = '/Users/pmm/Documents/xtal_learning/hydrocarbons/'

fnames= []
fnames.append(os.path.join(dir,'Data_ML_6.csv'))
fnames.append(os.path.join(dir,'Data_ML_7.csv'))
fnames.append(os.path.join(dir,'Data_ML_8.csv'))
fnames.append(os.path.join(dir,'Data_ML9.csv'))
fnames.append(os.path.join(dir,'Data_ML10.csv'))


manip = [("Name",func)]
usecols = ['Structures','Space group', 'Cell volume', 'Density', 
           'Total energy', 'van der Waals energy', 'Electrostatic energy']

rename ={'Structures':'Name',
         'Total energy':'Energy',
          'van der Waals energy': 'vdW Energy',
          'Electrostatic energy': 'Electrostatic Energy'}

df = csv_to_dataframe(fnames,manip=manip,rename=rename,index_col='Name',usecols=usecols)
df.to_pickle(os.path.join(dir,'Data_ML_set'))
print(df.head())
print(df.columns)