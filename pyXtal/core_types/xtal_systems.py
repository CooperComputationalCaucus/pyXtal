#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 16:58:19 2018

@author: alggroup
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:45:51 2018

@author: alggroup
"""

import numpy as np
from numpy import sqrt,cos,sin,degrees

__all__ = [
    'Xtal', 'No_lattice', 'zero_cubic', 'random_cubic'
    ]

class Xtal():
    """Builds the crystal as a list of atoms with associated
    dictionary for lattice paramters and 
    methods for conversions of basis etc.
    - atoms: list of atoms of specified data type
        - position
        - symbol
        - label
        - index
        - on
        - Uij
    - number_of_atoms: integer number of atoms
    - lat: dictionary of lattice parameters
    - bonding: dictionary of indicies for bonded members [a]=[b c d]
    - adjacency: adjacency matrix for bonded memebers. dim=(N_atoms,N_atoms)
    - isfract: logical dictating if the atomic coordinates are fractional
    - ac2f: matrix for left multiplication to go from cartesian to fractional coords
    - af2c: matrix for left multiplication to go from fractional to cartesian coords
    """
    def __init__(self,positions=None,lat=None,symbols=None,labels=None,indicies=None,isfract=True,bonding={},Uij=None):
        if positions is None:
            raise ValueError
        if symbols is None:
            symbols=["Si"]*len(positions)
        if labels is None:
            labels=[]
            for i,sym in enumerate(symbols):
                labels.append("{:s}{:d}".format(sym.strip(),i+1))
        if indicies is None:
            #indicies starting at 0 not 1.
            indicies=[]
            for i,sym in enumerate(symbols):
                indicies.append(i)
        self.number_of_atoms=len(positions)
        dt=atom_dt()
        self.atoms=np.zeros(self.number_of_atoms,dtype=dt)
        self.atoms['position']=positions
        self.atoms['label']=labels
        self.atoms['symbol']=symbols
        self.atoms['index']=indicies
        if Uij is None:
            self.atoms['Uij']=None
        else:
            self.atoms['Uij']=Uij
        self.atoms['on']=True
        self.lat=lat
        self.isfract=isfract
        self.af2c= np.zeros((3,3))
        if lat is None:
            raise Exception(("An instance of the Xtal class requires periodic boundary conditions "
                             "(i.e. lattice parameters) Consider using a different class."))
        try:
            v = sqrt(1 - cos(lat['alpha'])**2 - cos(lat['beta'])**2 - cos(lat['gamma'])**2 +
                 2*cos(lat['alpha'])*cos(lat['beta'])*cos(lat['gamma']))
            self.af2c[0,0] = self.lat['L_a']
            self.af2c[0,1] = self.lat['L_b']*cos(self.lat['gamma'])
            self.af2c[1,1] =self.lat['L_b']*sin(self.lat['gamma'])
            self.af2c[0,2] = self.lat['L_c']*cos(self.lat['beta'])
            self.af2c[1,2] = self.lat['L_c']*(cos(self.lat['alpha'])-cos(self.lat['beta'])*cos(self.lat['gamma']))/sin(self.lat['gamma'])
            self.af2c[2,2] = self.lat['L_c'] * v / sin(self.lat['gamma'])        
            self.ac2f = np.linalg.inv(self.af2c)
        except KeyError:
            raise KeyError("Key error in lattice parameters dictionary")
        except ValueError:
            raise ValueError("Value error in lattice parameters dictionary")
        #Presently implementing both a bonding diciotnary of indicies, and adjacency matrix to describe connectivity
        self.bonding=bonding
        self.adjacency=np.zeros((self.number_of_atoms,self.number_of_atoms))

    def __str__(self):
        return ("Crystal Lattice containing %i atoms.\n"
                "a = %f\nb = %f\nc = %f\nalpha = %f\nbeta = %f\ngamma = %f\n" %
                (len(self.atoms),self.lat['L_a'],self.lat['L_b'],self.lat['L_c'],
                 degrees(self.lat['alpha']),degrees(self.lat['beta']),degrees(self.lat['gamma'])))
        
    def recalc_conversions(self):
        lat= self.lat
        self.af2c= np.zeros((3,3))
        v = sqrt(1 - cos(lat['alpha'])**2 - cos(lat['beta'])**2 - cos(lat['gamma'])**2 +
             2*cos(lat['alpha'])*cos(lat['beta'])*cos(lat['gamma']))
        self.af2c[0,0] = self.lat['L_a']
        self.af2c[0,1] = self.lat['L_b']*cos(self.lat['gamma'])
        self.af2c[1,1] =self.lat['L_b']*sin(self.lat['gamma'])
        self.af2c[0,2] = self.lat['L_c']*cos(self.lat['beta'])
        self.af2c[1,2] = self.lat['L_c']*(cos(self.lat['alpha'])-cos(self.lat['beta'])*cos(self.lat['gamma']))/sin(self.lat['gamma'])
        self.af2c[2,2] = self.lat['L_c'] * v / sin(self.lat['gamma'])        
        self.ac2f = np.linalg.inv(self.af2c)
        
    def distance(self,atom1,atom2):
        """First works on fractional coordinates, then converts to cartesian distance
        Returns the magnitude of the cartesian distance, NOT the distance vector
        """
        self.cart2fract()
        dist = np.subtract(atom2['position'],atom1['position'])
        for i in range(dist.shape[0]):
            if abs(dist[i])>0.5:
                dist[i]=dist[i]-1*np.sign(dist[i])
        dist = np.dot(self.af2c,dist)
        return np.linalg.norm(dist)
    
    def configure_bonding(self,sym1,sym2,dmax):
        #Bonding uses the same indexing as Xtal, which starts at 0 by default. 
        for a in self.atoms:
            self.bonding[a['index']]=[]
            for b in self.atoms:
                if a['symbol'] == sym1 and b['symbol'] == sym2 and a['index'] != b['index']:
                    d = self.distance(a,b)
                    if d<=dmax:
                        self.bonding[a['index']].append(b['index'])
    
    def configure_adjacency(self):
        if self.bonding=={}:
            return
        for (a,bonded) in self.bonding.items():
            for b in bonded:
                self.adjacency[a,b]=1
                self.adjacency[b,a]=1
    
    def fract2cart(self):
        if not self.isfract:
            return
        
        for atom in self.atoms:
            atom['position'] = np.dot(self.af2c,atom['position'])            
        self.isfract = False
        
    def cart2fract(self):
        if self.isfract:
            return
        
        for atom in self.atoms:
            atom['position'] = np.dot(self.ac2f,atom['position']) 
        self.isfract  = True
    
class No_lattice(): 
    """Builds system box with no periodic boundary conditions.
    Always considered to be in cartesian coordinates, and never fractional coordinates.
    - atoms: list of atoms of specified data type
    """
    def __init__(self,positions=None,symbols=None,labels=None,indicies=None,bonding={}):
        if positions is None:
            raise ValueError
        if symbols is None:
            symbols=["Si"]*len(positions)
        if labels is None:
            labels=[]
            for i,sym in enumerate(symbols):
                labels.append("{:2s}{:03d}".format(sym,i+1))
        if indicies is None:
            indicies=[]
            for i,sym in enumerate(symbols):
                indicies.append(i+1)
        self.number_of_atoms=len(positions)
        dt=atom_dt()
        self.atoms=np.zeros(self.number_of_atoms,dtype=dt)
        self.atoms['position']=positions
        self.atoms['label']=labels
        self.atoms['symbol']=symbols
        self.atoms['index']=indicies
        self.atoms['on']=True
        self.bonding=bonding

def atom_dt():
    """Builds the data type for the atoms, this consists of:

    - position
    - label
    - symbol
    - index
    - on (Logical)
    """
    return np.dtype([('label', np.str_, 16),('position',np.float64,(3,)),('symbol',np.str_,2),
                     ('index',np.int),('on',bool),('Uij',np.float64,(3,3))])
    
def zero_cubic(N_atoms):
    '''Returns a simple cubic system with all silicon atoms at position [0,0,0]. 
    Mainly a function for testing, but could be used for system initialization. 
    '''
    dt=np.dtype((np.float64,(3,)))
    pos=np.zeros(N_atoms,dt)
    cell_data={}
    cell_data['L_a'] =  1.0
    cell_data['L_b'] = 1.0
    cell_data['L_c'] = 1.0
    cell_data['alpha'] = np.radians(90.)
    cell_data['beta'] = np.radians(90.0)
    cell_data['gamma'] = np.radians(90.0)
    symbols=["Si"]*N_atoms
    xtal=Xtal(positions=pos,lat=cell_data,symbols=symbols)
    
    return xtal

def random_cubic(N_atoms):
    '''Returns a simple cubic system with all silicon atoms at random positions. 
    Mainly a function for testing, but could be used for system initialization. 
    '''
    pos=np.random.rand(N_atoms,3)
    cell_data={}
    cell_data['L_a'] =  1.0
    cell_data['L_b'] = 1.0
    cell_data['L_c'] = 1.0
    cell_data['alpha'] = np.radians(90.)
    cell_data['beta'] = np.radians(90.0)
    cell_data['gamma'] = np.radians(90.0)
    symbols=["Si"]*N_atoms
    xtal=Xtal(positions=pos,lat=cell_data,symbols=symbols)
    
    return xtal

