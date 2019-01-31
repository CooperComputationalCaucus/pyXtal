#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Requires spglib and ASE to be installed. 



@author: PMM
"""

import ase
import ase.io
import spglib

def findsym(fname,prec=None,verbose=True):
    '''
    Reads cif file using ASE parser and extracts symmetry using spglib
    over different levels of tolerence
    '''
    
    if prec:
        assert(type(prec)==type(1e-10))
        precs=[prec]
    else:
        t = [1, 2, 3, 5, 7]
        precs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6] + [i * j for i in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1] for j in t]
    
    old = ''
    if verbose: print('# Tolerance\tSpace group')
    
    cell = ase.io.read(fname)
    for prec in precs:
        s = spglib.get_spacegroup(cell, symprec=prec)
        if s != old:
            if verbose: print(f'{prec}\t\t{s}')
            old = s

    sg, sgid = old.split(' (')
    sgid = sgid.replace(')', '')

    return sg,sgid

def sort_cifs(dir,precision=None,copy=True):
    '''
    Takes all of the cifs in a directory and places them in subfolders according to their 
    spacegroup number.  
    
    If no precision 
    
    Inputs
    ======================
    dir: directory of cifs
    precision : precision of symmetry calculation. If none is given the highest symmetry est is given
    copy: True to create copies of files, False to move files
    '''
    import os
    from shutil import copyfile, move
    
    file_list = os.listdir(dir)

    for file in file_list:
        if os.path.splitext(file)[1] == '.cif':
            inpath = os.path.join(dir,file)
            print(file)
            num = findsym(inpath,verbose=False)[1]+'/'
            symdir = os.path.join(dir,num)

            if not os.path.isdir(symdir):
                os.mkdir(symdir)
            outpath = os.path.join(symdir,file)

            if copy:
                copyfile(inpath,outpath)
            else:
                move(inpath,outpath)
    
if __name__ == '__main__':
    sort_cifs('/Users/pmm/Documents/xtal_learning/T2_findsym/cifs',copy=False)
    
