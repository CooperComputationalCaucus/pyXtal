#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:46:26 2018

@author: alggroup
"""
import numpy as np
import sys
import pyXtal as pxt

#dt=np.dtype((np.float64,(3,)))
#pos=np.zeros(5,dt)
#cell_data={}
#cell_data['L_a'] =  1.
#cell_data['L_b'] = 1.0
#cell_data['L_c'] = 1.0
#cell_data['alpha'] = radians(90.)
#cell_data['beta'] = radians(90.0)
#cell_data['gamma'] = radians(90.0)
#try:
#    newXtal=pxt.Xtal(positions=pos,lat=cell_data)
#    print(newXtal.atoms)
#except:
#    print(sys.exc_info()[0])
#    print(sys.exc_info()[1])
#try:
#    newXtal=pxt.random_cubic(5)
#    print(newXtal.atoms)    
#except:
#    print(sys.exc_info()[0])
#    print(sys.exc_info()[1])

pxt.io_types.io_tests()

#if isinstance(newXtal,Xtal):
#    print('yes')    
#    print(newXtal)
#else:
#    print('no')       
#        
