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

def percent_energy_difference(df):
    '''
    Returns the percent energy difference from the global minimum as a df column
    '''
    min = np.min(df['Energy'])
    C = (min - df['Energy'])/min*100
    return C

class feature_scaling():
    '''
    Class of feature scaling that requires specific feature in constructor.
    
    Mainly for working with dataframes. 
    '''
    def __init__(self,feature):
        assert type(feature)==type(" ")
        self.feature = feature
         
    def mean_scaling(self,df):
        '''
        Returns new series of mean scaled feature
        '''
        min_f = np.min(df[self.feature])
        max_f = np.max(df[self.feature])
        mean_f = np.mean(df[self.feature])
        return (df[self.feature]- mean_f)/(max_f-min_f)
   
    def min_max_scaling(self,df):
        '''
        Returns new series of min-max scaled feature
        '''
        min_f = np.min(df[self.feature])
        max_f = np.max(df[self.feature])
        return (df[self.feature]- min_f)/(max_f-min_f)
    
    
def csv_to_dataframe(fnames,manip=[],rename={},index_col=None,usecols=None,derived_cols=None):
    '''
    Imports csv file and outputs dataframe.
    If a list of file names are give, the resultant dataframe is stacked. 
    
    Inputs
    ==================    
    fnames : CSV path or list of csv paths
    manip: list of tuples of inplace colummn labels/manipulation function pairs (functions act on column items)
    rename : dictionary of renaming labels from keys to values
    index_col:  NEW label to index by. Keeps default indexing for None. 
    usecols: columns to use from the import CSV file
    derived_cols: list of tuples of new colummn labels/manipulation function pairs (functions act of dataframe)
    
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
        # Run column derivations
        for label, func in derived_cols:
            _df[label] = func(_df)
        #Set index
        if index_col: _df.set_index(index_col,inplace=True)
        frames.append(_df)
    
    df = pd.concat(frames)
    return df

def parse_filenames(dir,keys=[],sep='_'):
    '''
    Parses filenames in a directory, and separates them into keys as presented, 
    and outputs a dataframe. If '_' is given as a key, the data is discarded.
    If not enough keys are given, only the first components are assigned. 
    
    Inputs
    ==================    
    dir : directory of interest
    sep : string seperator
    keys : list of keys for data splitting. 
    
    Returns
    ==================
    df : complete dataframe
    '''
    data = {}
    for key in keys:
        if key == '_': continue
        data[key] = []
    data["Name"] = []
    fnames = os.listdir(dir)
    for fname in fnames:
        base = os.path.splitext(fname)[0]
        data["Name"].append(base)
        for idx,col in enumerate(base.split(sep)[:len(keys)]):
            if keys[idx] == '_': continue
            data[keys[idx]].append(col)
                
    df = pd.DataFrame.from_dict(data)      
    return df