import numpy as np
import sys
sys.path.append('../')
import pyXtal as pxt
xtal=pxt.read_cif('./forzac.cif')
pxt.write_cif('./forzac_new.cif',xtal,True)
print(xtal)