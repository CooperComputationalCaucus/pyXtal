'''
Script to pull in set of cif files and make a single dataframe .

@author: pmm
'''
import numpy as np
import os
import sys
sys.path.append('../')
import pyXtal as pxt
from pyXtal.csp_utils.dataset_io import parse_filenames
import pandas as pd

dirs = ['/Users/pmm/Documents/xtal_learning/triptycene/cifs/T2',
        '/Users/pmm/Documents/xtal_learning/triptycene/cifs/ring3',
        '/Users/pmm/Documents/xtal_learning/triptycene/cifs/ring32',
        '/Users/pmm/Documents/xtal_learning/triptycene/cifs/ring34',
        '/Users/pmm/Documents/xtal_learning/triptycene/cifs/ring39'
        ]
frames=[]
for dir in dirs:
    _df = parse_filenames(dir, keys=['_','Energy','Density'])
    _df.set_index("Name",inplace=True)
    frames.append(_df)
    
df = pd.concat(frames)
df.to_pickle('/Users/pmm/Documents/xtal_learning/triptycene/cifs/data_triptycene_set')