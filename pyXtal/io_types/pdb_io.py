#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 17:30:26 2018

@author: Phillip M. Maffettone
"""
from ..core_types import *
import os
import numpy as np
__all__ = [
    'read_pdb', 'write_pdb'
    ]

def read_pdb(filename, pbc=False):
    '''
    Routine for reading .pdb files.
    Specifically works for atom chains with connectivity. 
    This does not interpret a vast majority of the information in a pdb file,
    but is suffiecient for the necessary implementation in pymol. 
    Places pdb file in new style dictionary and returns an Xtal class
    Presently assumes the indexes run from 1 to num_of_atoms. 
    ###############
    Parameters:
        filename - complete location of pdb file
    ###############
    '''
    with open(filename, 'r') as f:
        line_num = 0
        x=[]
        y=[]
        z=[]
        idxs=[]
        syms=[]
        bonding={}
        for line in f:
            line_num += 1
            cols = line.split()
            if cols[0] ==  'HEADER':
                continue
            elif cols[0] == 'ATOM':
                idxs.append(int(cols[1]))
                syms.append(cols[2])
                x.append(float(cols[6]))
                y.append(float(cols[7]))
                z.append(float(cols[8]))
            elif cols[0] == 'CONECT':
                bonding[int(cols[1])-1]=[]
                for idx in cols[2:]:
                    bonding[int(cols[1])-1].append(int(idx)-1)
    coords=np.transpose(np.array([x,y,z]))

    if pbc:
        cell_data={}
        x_min=+1e6
        x_max=-1e6
        y_min=+1e6
        y_max=-1e6
        z_min=+1e6
        z_max=-1e6
        for i in coords[:,0]:
            if i<x_min:
                x_min=i
            elif i>x_max:
                x_max=i
        for i in coords[:,1]:
            if i<y_min:
                y_min=i
            elif i>y_max:
                y_max=i
        for i in coords[:,2]:
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
        xtal=Xtal(positions=coords,lat=cell_data,symbols=syms,isfract=False,bonding=bonding)
    else:
        xtal = No_lattice(positions=coords,symbols=syms,bonding=bonding)
                
    return xtal
def write_pdb(filename,xtal,connective=False,redundant_boundary=False,box=False,verbose=False):
    '''
    Routine for writing PDB file to include bonding
    '''
    xtal.frac2cart()
    if connective:
        write_pdb_connective(filename,xtal,redundant_boundary=redundant_boundary,verbose=verbose)
    else:
        write_pdb_truncated(filename,xtal,redundant_boundary=redundant_boundary,verbose=verbose)
    if unit_cell:
        write_pdb_box(filename,xtal,verbose)
        
    
def write_pdb_truncated(filename,xtal,redundant_boundary=False,verbose=False):
    '''
    Routine for writing PDB file to include bonding. Does not include bonds across PBC
    '''
    atoms = xtal.atoms
    fmtstr= "{:6s}{:5d}{:1s}{:4s}{:1s}{:3s}{:1s}{:1s}{:4d}{:1s}{:3s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:10s}{:>2s}{:2s}\n"
    with open(filename,'w') as f:
        f.write('HEADER    pyXtal generation\n')
        for i,p in enumerate(atoms['position']):
            f.write(fmtstr.format('ATOM',atoms['index'][i]+1,'',atoms['symbol'][i],'','MOL','','H',0,'','',p[0], 
                                  p[1],p[2],1.0,0.0,'',atoms['symbol'][i],'0'))
        xtal.cart2fract()
        for (a,bonded) in xtal.bonding.items():
            f.write('{:6s}{:5d}'.format('CONECT',a+1))
            for b in bonded:
                if (all(x < 0.5 for x in abs(np.subtract(atoms['position'][a],atoms['position'][b])))):
                    f.write('{:5d}'.format(b+1))
            f.write('\n')               
        f.write('TER')
    if verbose:
        xtal.fract2cart()
        atoms = xtal.atoms
        fmtstr= "{:6s}{:5d}{:1s}{:4s}{:1s}{:3s}{:1s}{:1s}{:4d}{:1s}{:3s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:10s}{:>2s}{:2s}"
        print('HEADER    pyXtal generation')
        for i,p in enumerate(atoms['position']):
            print(fmtstr.format('ATOM',atoms['index'][i]+1,'',atoms['symbol'][i],'','MOL','','H',0,'','',p[0], 
                                  p[1],p[2],1.0,0.0,'',atoms['symbol'][i],'0'))
        xtal.cart2fract()
        for (a,bonded) in xtal.bonding.items():
            print('{:6s}{:5d}'.format('CONECT',a+1), end="")
            for b in bonded:
                if (all(x < 0.5 for x in abs(np.subtract(atoms['position'][a],atoms['position'][b])))):
                    print('{:5d}'.format(b+1),end="")
            print("")
        print("TER")
        
def write_pdb_connective(filename,xtal,redundant_boundary=False,verbose=False):
    '''
    Routine for writing PDB file to include all bonding
    '''
    atoms = xtal.atoms
    fmtstr= "{:6s}{:5d}{:1s}{:4s}{:1s}{:3s}{:1s}{:1s}{:4d}{:1s}{:3s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:10s}{:>2s}{:2s}\n"
    with open(filename,'w') as f:
        f.write('HEADER    pyXtal generation\n')
        for i,p in enumerate(atoms['position']):
            f.write(fmtstr.format('ATOM',atoms['index'][i]+1,'',atoms['symbol'][i],'','MOL','','H',0,'','',p[0], 
                                  p[1],p[2],1.0,0.0,'',atoms['symbol'][i],'0'))
        xtal.cart2fract()
        for (a,bonded) in xtal.bonding.items():
            f.write('{:6s}{:5d}'.format('CONECT',a+1))
            for b in bonded:
                f.write('{:5d}'.format(b+1))
            f.write('\n')               
        f.write('TER')
    if verbose:
        xtal.fract2cart()
        atoms = xtal.atoms
        fmtstr= "{:6s}{:5d}{:1s}{:4s}{:1s}{:3s}{:1s}{:1s}{:4d}{:1s}{:3s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:10s}{:>2s}{:2s}"
        print('HEADER    pyXtal generation')
        for i,p in enumerate(atoms['position']):
            print(fmtstr.format('ATOM',atoms['index'][i]+1,'',atoms['symbol'][i],'','MOL','','H',0,'','',p[0], 
                                  p[1],p[2],1.0,0.0,'',atoms['symbol'][i],'0'))
        xtal.cart2fract()
        for (a,bonded) in xtal.bonding.items():
            print('{:6s}{:5d}'.format('CONECT',a+1), end="")
            for b in bonded:
                print('{:5d}'.format(b+1),end="")
            print("")
        print("TER")
        
def write_pdb_box(filename,xtal,verbose=False):
    filename, file_extension = os.path.splitext(filename)
    filename=filename+'_unitcell'+file_extension
    coords = np.array([[0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0],
                        [1.0, 1.0, 0.0],
                        [1.0, 0.0, 1.0],
                        [0.0, 1.0, 1.0],
                        [1.0, 1.0, 1.0],])
    syms=['Ar','Ar','Ar','Ar','Ar','Ar','Ar']
    box=Xtal(positions=coords,lat=xtal.lat,symbols=syms)
    box.configure_bonding('Ar','Ar',dmax=(box.lat['L_a']+.01))
    write_pdb_truncated(filename,box,redundant_boundary=False,verbose=verbose)
