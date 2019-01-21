#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 17:10:10 2018

@author: alggroup
"""

import numpy as np
import sys
sys.path.append('../')
import pyXtal as pxt

xtal=pxt.read_cif('./testcif2.cif')
xtal.configure_bonding('Si','Si', 4.5)
print(xtal.bonding)

#pxt.write_cif('./testout.cif',xtal,verbose=True)
#xtal=pxt.read_xyz("./testxyz.xyz",pbc=True)
#xtal=pxt.read_pdb("./testpdb.pdb",pbc=True)
#pxt.write_xyz('./testout.xyz',xtal,verbose=True)
#pxt.write_cif('./testout.cif',xtal,verbose=True)
pxt.write_pdb('./testout.pdb',xtal,verbose=True)
pxt.write_pdb_connective('./testout.pdb',xtal,verbose=True)