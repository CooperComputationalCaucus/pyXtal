"""
Module to handle res datafiles from mercury and CSP

@author: PMM
"""

from ..core_types import *
from numpy import radians,degrees
import numpy as np
import os

__all__ = [
    'read_res' 
    ]

def read_res(filename):
    '''
    Routine for reading .RES files.
    Built off of a presumed format from a set of files sent from CSP datasets.
    
    Inputs
    ==============
    filename: path to .RES file to read
    
    TODO:
     - Actually write this
    '''
    