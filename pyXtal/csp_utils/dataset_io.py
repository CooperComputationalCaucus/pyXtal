#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 17:30:26 2018

Some brief functions to handle csp datasets

@author: PMM
"""

import numpy as np
import os
import pandas as pd

def replace_spaces(string):
    string = '_'.join(string.split(' '))
    return string

def csv_to_dataframe(fnames,manip=[],rename={},index_col=None,usecols=None):
    '''
    Imports csv file and outputs dataframe.
    If a list of file names are give, the resultant dataframe is stacked. 
    
    Inputs
    ==================    
    fnames : CSV path or list of csv paths
    manip: list of tuples of new column labels/manipulation function pairs
    rename : dictionary of renaming labels from keys to values
    index_col:  NEW label to index by. Keeps default indexing for None. 
    usecols: columns to use from the import CSV file
    
    Returns
    ==================
    df : complete dataframe
    '''
    
    if type(fnames) != type([]):
        fnames = [fnames]
    
    frames = []
    for fname in fnames:
        # Read csv
        _df = pd.read_csv(fname,usecols=usecols)
        # Change necessary column labels
        _df.rename(index=str,columns=rename,inplace=True)
        # Run manipulaitons
        for label,func in manip:
            _df[label] = _df[label].apply(func)
        #Set index
        if index_col: _df.set_index(index_col,inplace=True)
        frames.append(_df)
    
    df = pd.concat(frames)
    return df