#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 17:30:26 2018

@author: alggroup
"""

from ..core_types import *
from numpy import radians,degrees
import numpy as np

__all__ = [
    'read_xyz', 'write_xyz', 'xyz_io_test', 
    ]


def read_xyz(filename,symbol='Si',pbc=False):
    '''
    Reads a .xyz file where the first line contains the number of atoms
    and the subsequent lines contain 3-vectors for cartesian coordinates
    The I/O will handle atomic symbols leading the coordinates, but does not handle 
    any comment lines. 
    ###############
    Parameters:
    filename - complete location of .xyz file
    symbol - atomic symbol to give every coordinate 
    pbc - logical to dictate addition of artifical periodic boundary conditions
    ###############
    '''
    with open(filename,'r') as f:
        line_num = 0
        x=[]
        y=[]
        z=[]
        symbols=[]
        for line in f:
            line_num +=1
            cols = line.split()
            if len(cols)==0: continue
            if line_num==1: continue
            try:
                float(cols[0])
                x.append(float(cols[0].split('(',1)[0]))
                y.append(float(cols[1].split('(',1)[0]))
                z.append(float(cols[2].split('(',1)[0]))
            except(ValueError):
                symbols.append(cols[0])
                x.append(float(cols[1].split('(',1)[0]))
                y.append(float(cols[2].split('(',1)[0]))
                z.append(float(cols[3].split('(',1)[0]))
    if not symbols: 
        symbols=[symbol]*len(x)
    elif len(x) != len(symbols):
        raise InputError("The read number of atomic symbols and positions does not match. Check input file.\n")

    coords=np.transpose(np.array([x,y,z]))
    
    if pbc:
        cell_data={}
        x_min=+1e6
        x_max=-1e6
        y_min=+1e6
        y_max=-1e6
        z_min=+1e6
        z_max=-1e6
        for i in x:
            if i<x_min:
                x_min=i
            elif i>x_max:
                x_max=i
        for i in y:
            if i<y_min:
                y_min=i
            elif i>y_max:
                y_max=i
        for i in z:
            if i<z_min:
                z_min=i
            elif i>z_max:
                z_max=i
        cell_data['L_a'] = x_max-x_min
        cell_data['L_b']= y_max-y_min
        cell_data['L_c']= z_max-z_min
        cell_data['alpha'] = radians(90.0)
        cell_data['beta'] = radians(90.0)
        cell_data['gamma'] = radians(90.0)
        xtal = Xtal(positions=coords,lat=cell_data,symbols=symbols,isfract=False)
    else:
        xtal = No_lattice(positions=coords,symbols=symbols)
    return xtal
        

def write_xyz(filename, xtal, verbose=False):
    '''
    Routine for writing an .xyz file
    Can handle a crystal with or without pbc
    ###############
    Parameters:
    filename - complete location of .xyz file
    xtal - system of type Xtal or No_lattice
    verbose - logical for printing .xyz representation to buffer
    ###############
    '''
    if isinstance(xtal,Xtal):
        xtal.fract2cart()
    with open(filename,'w') as f:
        f.write('{:6d}  {:2s}'.format(len(xtal.atoms['symbol']),xtal.atoms['symbol'][0]))
        for i in range(len(xtal.atoms['symbol'])):
            f.write("{:12.8f}{:12.8f}{:12.8f}".format(
                xtal.atoms['position'][i][0],xtal.atoms['position'][i][1],xtal.atoms['position'][i][2]))
            if verbose:
                print("{:12.8f}{:12.8f}{:12.8f}".format(
                    xtal.atoms['position'][i][0],xtal.atoms['position'][i][1],xtal.atoms['position'][i][2]))

def xyz_io_test():
    ''' 
    Quick and dirty test to make sure module is loaded
    '''
    return True
