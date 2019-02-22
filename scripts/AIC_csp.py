'''
Script to pull in set of csv files and make on single dataframe .

@author: pmm
'''
import numpy as np
import os
import sys
sys.path.append('../')
import pyXtal as pxt
from pyXtal.csp_utils.dataset_io import csv_to_dataframe, feature_scaling

func = pxt.csp_utils.dataset_io.replace_spaces
func2 = pxt.csp_utils.dataset_io.percent_energy_difference
f_scaling = feature_scaling("Energy")

dir = '/Users/pmm/Documents/xtal_learning/hydrocarbons/'

fnames= []
fnames.append(os.path.join(dir,'Data_ML_1.csv'))
fnames.append(os.path.join(dir,'Data_ML_2.csv'))
fnames.append(os.path.join(dir,'Data_ML_3.csv'))
fnames.append(os.path.join(dir,'Data_ML_4.csv'))
fnames.append(os.path.join(dir,'Data_ML_5.csv'))
fnames.append(os.path.join(dir,'Data_ML_6.csv'))
fnames.append(os.path.join(dir,'Data_ML_7.csv'))
fnames.append(os.path.join(dir,'Data_ML_8.csv'))
fnames.append(os.path.join(dir,'Data_ML_9.csv'))
fnames.append(os.path.join(dir,'Data_ML_10.csv'))
fnames.append(os.path.join(dir,'Data_ML_35.csv'))
fnames.append(os.path.join(dir,'Data_ML_49.csv'))
fnames.append(os.path.join(dir,'Data_ML_50.csv'))
fnames.append(os.path.join(dir,'Data_ML_54.csv'))




manip = [("Name",func)]
derived_cols = [("PercentED",func2),
                ("min_max_scaled_E", f_scaling.min_max_scaling),
                ("mean_scaled_E", f_scaling.mean_scaling)]
usecols = ['Structures','Space group', 'Cell volume', 'Density', 
           'Total energy', 'van der Waals energy', 'Electrostatic energy']

rename ={'Structures':'Name',
         'Total energy':'Energy',
          'van der Waals energy': 'vdW Energy',
          'Electrostatic energy': 'Electrostatic Energy'}

df = csv_to_dataframe(fnames,manip=manip,rename=rename,index_col='Name',usecols=usecols,derived_cols=derived_cols)

df.to_pickle(os.path.join(dir,'Data_ML_set'))
df.to_csv(os.path.join(dir,'Data_ML_set.csv'))
print(df.columns)
print(df.index)
