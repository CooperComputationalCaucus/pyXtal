#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Script to split CSP cifs from Chengxi and log density and energy data into dataframe pickle

@author: pmm
'''
import numpy as np
import os
import sys
sys.path.append('../')
import pyXtal as pxt

dir = '/Users/pmm/Documents/xtal_learning/triptycene_chengxi/'
fname = 'totalT2fixed.cif'
dir2 = os.path.join(dir,os.path.splitext(fname)[0])
os.mkdir(dir2)

pxt.split_cif(os.path.join(dir,fname),dir=dir2)