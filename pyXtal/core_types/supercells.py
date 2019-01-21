#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:32:40 2018

@author: Phillip M. Maffettone
"""

from ..core_types import *
import numpy as np
__all__ = [
    'build_super'
    ]

def build_super(xtal,Nx=1,Ny=1,Nz=1,bonded=False):
    '''
    Builds an Nx x Ny x Nz super cell from a 
    list of atoms in fractional coordinates.
    Rewrites atom indexing.
    '''
    xtal.cart2fract()
    if bonded:
        build_super_bonded(xtal,Nx=1,Ny=1,Nz=1)
    else:
        build_super_unbonded(xtal,Nx=1,Ny=1,Nz=1)

def build_super_unbonded(xtal,Nx=1,Ny=1,Nz=1):
    pass

def build_super_bonded(xtal,Nx=1,Ny=1,Nz=1):
    pass