#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 17:30:26 2018

@author: PMM
"""

from ..core_types import *
from numpy import radians,degrees
import numpy as np
import os
import numexpr

__all__ = [
    'read_cif', 'write_cif','split_cif', 
    ]

def read_cif(filename):
    '''
    Routine for reading .CIF files.
    Structured as a first attempt with i/o in python, using many comments
    Places cif file in new style dictionary and returns an Xtal class
    Presently requires plain text descriptions of symetry ops, and does not interpret space groups
    ###############
    Parameters:
        filename - complete location of cif file
    ###############
    '''
    data = {}
    cell_data = {}
    
    # Open the .CIF file
    with open(filename,'r') as f:
        reading_sym_ops = False
        reading_atom_sites = False
        reading_atom_aniso = False
        line_num = 0

        # Read lines one by one
        for line in f:
            line_num += 1

            # Split into columns
            cols = line.split()
            if (len(cols) == 0): continue

            # ID the keywords
            # Keywords associated with cell parameters
            if (cols[0] == '_cell_length_a'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_a'] = float(cols[1])
            elif (cols[0] == '_cell_length_b'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_b'] = float(cols[1])
            elif (cols[0] == '_cell_length_c'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_c'] = float(cols[1])
            elif (cols[0] == '_cell_angle_alpha'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_alpha'] = float(cols[1])
            elif (cols[0] == '_cell_angle_beta'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_beta'] = float(cols[1])
            elif (cols[0] == '_cell_angle_gamma'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_gamma'] = float(cols[1])
            elif (cols[0] == '_cell_volume'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_volume'] = float(cols[1])

            # Keywords associated with symmetry operations
            elif (cols[0] == '_symmetry_equiv_pos_as_xyz'):
                reading_sym_ops = True
                data['_symmetry_equiv_pos_as_xyz'] = []
            elif (reading_sym_ops):
                # Add the operation if the string is between single quotes.
                # Otherwise it's a sign we are done with the list.
                if (cols[0][0] == '\''  and  cols[0][-1] == '\''):
                    data['_symmetry_equiv_pos_as_xyz'].append(cols[0][1:-1])
                elif(cols[0].isdigit()):
                    data['_symmetry_equiv_pos_as_xyz'].append(cols[1])
                else:
                    reading_sym_ops = False

            # Keywords associated with atom sites
            elif (cols[0] == '_atom_site_label'):
                index_label = line_num
                data['_atom_site_label'] = []
                data['_atom_site_fract_x'] = []
                data['_atom_site_fract_y'] = []
                data['_atom_site_fract_z'] = []
                reading_atom_sites = True

            # Keep track of where the other labels are (order is important).
            elif (cols[0] == '_atom_site_fract_x'):
                index_x = line_num - index_label
            elif (cols[0] == '_atom_site_fract_y'):
                index_y = line_num - index_label
            elif (cols[0] == '_atom_site_fract_z'):
                index_z = line_num - index_label
            elif (cols[0] == '_atom_site_type_symbol'):
                index_sym = line_num - index_label
                data['_atom_site_type_symbol'] = []

            # If we are currently reading the atom sites...
            elif (reading_atom_sites):
                # Read the actual data if we have 4 columns or more of data.
                if (len(cols) >= 4):
                    data['_atom_site_label'].append(cols[0])
                    data['_atom_site_fract_x'].append(float(cols[index_x].split('(',1)[0]))
                    data['_atom_site_fract_y'].append(float(cols[index_y].split('(',1)[0]))
                    data['_atom_site_fract_z'].append(float(cols[index_z].split('(',1)[0]))
                    if ('_atom_site_type_symbol' in data): data['_atom_site_type_symbol'].append(cols[index_sym])

                # Stop reading atom sites if we found a line with fewer
                # columns, and which does not start with '_atom_site_'.
                elif (len(cols[0]) < 11  or  cols[0][:11] != '_atom_site_'):
                    reading_atom_sites = False
            
            #Keywords associated with aniso. Presently only works for isoptropic case with symmetry. Anisotropic without
            elif(cols[0] == '_atom_site_aniso_label'):
                index_label = line_num
                data['_atom_site_aniso_label'] = []
                data['_atom_site_aniso_U_11'] = []
                data['_atom_site_aniso_U_22'] = []
                data['_atom_site_aniso_U_33'] = []
                data['_atom_site_aniso_U_23'] = []
                data['_atom_site_aniso_U_13'] = []
                data['_atom_site_aniso_U_12'] = []
                reading_atom_aniso = True

            # Keep track of where the other labels are (order is important).
            elif(cols[0] == '_atom_site_aniso_U_11'):
                index_11 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_22'):
                index_22 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_33'):
                index_33 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_23'):
                index_23 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_13'):
                index_13 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_12'):
                index_12 = line_num-index_label
            elif(reading_atom_aniso):
                # Read the actual data if we have 4 columns or more of data.
                if (len(cols) >= 4):
                    data['_atom_site_aniso_label'].append(cols[0])
                    data['_atom_site_aniso_U_11'].append(float(cols[index_11]))
                    data['_atom_site_aniso_U_22'].append(float(cols[index_22]))
                    data['_atom_site_aniso_U_33'].append(float(cols[index_33]))
                    data['_atom_site_aniso_U_23'].append(float(cols[index_23]))
                    data['_atom_site_aniso_U_13'].append(float(cols[index_13]))
                    data['_atom_site_aniso_U_12'].append(float(cols[index_12]))

                # Stop reading atom sites if we found a line with fewer
                # columns, and which does not start with '_atom_site_'.
                elif (len(cols[0]) < 11  or  cols[0][:11] != '_atom_site_'):
                    reading_atom_aniso = False
    
    #Creating a cleaner dictionary to return with a complete atom list
    cell_data['L_a'] =  float(data['_cell_length_a'])
    cell_data['L_b'] = float(data['_cell_length_b'])
    cell_data['L_c'] = float(data['_cell_length_c'])
    cell_data['alpha'] = radians(float(data['_cell_angle_alpha']))
    cell_data['beta'] = radians(float(data['_cell_angle_beta']))
    cell_data['gamma'] = radians(float(data['_cell_angle_gamma']))
    if ('_cell_volume' in data): cell_data['volume'] = float(data['_cell_volume'])

    #Creating atom list
    syms = data['_atom_site_type_symbol']
    labels = data['_atom_site_label']
    coords = np.transpose(np.array([data['_atom_site_fract_x'],data['_atom_site_fract_y'],
                       data['_atom_site_fract_z']]))
    
    if ('_atom_site_aniso_label' in data):
        Uij = np.transpose(np.array([
                [data['_atom_site_aniso_U_11'],data['_atom_site_aniso_U_12'],data['_atom_site_aniso_U_13']],
                [data['_atom_site_aniso_U_12'],data['_atom_site_aniso_U_22'],data['_atom_site_aniso_U_23']],
                [data['_atom_site_aniso_U_13'],data['_atom_site_aniso_U_23'],data['_atom_site_aniso_U_33']]
                        ]))
    else:
        Uij=None

    # Fixing with float division and applying symmetry operations
    if ('_symmetry_equiv_pos_as_xyz' in data):
        ops = data['_symmetry_equiv_pos_as_xyz']
        for i in range(len(ops)):
            ops[i] = ops[i].replace("/","./")

        # Two atoms are on top of each other if they are less than "eps" away.
        eps = 0.01  # in Angstrom
        imax = len(syms)
        new_coords=[]
        new_syms=[]
        new_labels=[]
        new_U=[]
        for i in range(imax):
            for op in ops:
                x,y,z = coords[i,:]
                # Text evaluation of symmetry operation (numxpr faield with list comprehensions)
                xn = numexpr.evaluate(op.split(',')[0]).item()
                yn = numexpr.evaluate(op.split(',')[1]).item()
                zn = numexpr.evaluate(op.split(',')[2]).item()
                
                # Forcing into the unit cell
                xn = (xn + 10.0) % 1.0
                yn = (yn + 10.0) % 1.0
                zn = (zn + 10.0) % 1.0
                # Adding only unique new atoms
                new_atom = True
                for c in range(imax):
                    if (abs(coords[c,0]-xn)<eps and abs(coords[c,1]-yn)<eps and abs(coords[c,2]-zn)<eps):
                        new_atom = False
                        break
                for c in range(len(new_coords)):
                    if (abs(new_coords[c][0]-xn)<eps and abs(new_coords[c][1]-yn)<eps and abs(new_coords[c][2]-zn)<eps):
                        new_atom = False
                        break
                if new_atom:
                    new_coords.append(np.array([xn,yn,zn]))
                    new_syms.append(syms[i])
                    new_labels.append(labels[i])
                    if Uij: new_U.append(Uij[i,:])

    #Sorting the atom list alphabetically
    if (new_coords):
        coords=np.vstack((coords, new_coords))
        syms.extend(new_syms)
        labels.extend(new_labels)
        if Uij: Uij=np.vstack((Uij,new_U))
    xtal = Xtal(positions=coords,lat=cell_data,symbols=syms,labels=labels,Uij=Uij)
    return xtal
def write_cif(filename,xtal,verbose=False):
    '''
    Routine for writing cif file in P1 
    from lattice parameters and atom list.
    Works directly from the returned objects of read_cif()
    ###############
    Parameters:
    filename - complete location of cif file
    xtal - crystal of type Xtal
    verbose - logical for printing .cif representation to buffer
    ###############
    '''
    xtal.cart2fract()
    atoms = xtal.atoms
    lat = xtal.lat
    with open(filename,'w') as f:
        f.write('data_cif\n\n')
        f.write('{:32s}{:7.4f}(0)\n'.format('_cell_length_a',lat['L_a']))
        f.write('{:32s}{:7.4f}(0)\n'.format('_cell_length_b',lat['L_b']))
        f.write('{:32s}{:7.4f}(0)\n'.format('_cell_length_c',lat['L_c']))
        f.write('{:32s}{:7.4f}(0)\n'.format('_cell_angle_alpha',degrees(lat['alpha'])))
        f.write('{:32s}{:7.4f}(0)\n'.format('_cell_angle_beta',degrees(lat['beta'])))
        f.write('{:32s}{:7.4f}(0)\n\n'.format('_cell_angle_gamma',degrees(lat['gamma'])))
        f.write("{:35s} 'P1'\n".format('_symmetry_space_group_name_H-M'))
        f.write('{:35s} 1\n'.format('_symmetry_Int_Tables_number'))
        f.write('{:35s} triclinic\n\n'.format('_symmetry_cell_setting'))
        f.write('loop_\n')
        f.write('_atom_site_label\n')
        f.write('_atom_site_type_symbol\n')
        f.write('_atom_site_fract_x\n')
        f.write('_atom_site_fract_y\n')
        f.write('_atom_site_fract_z\n')
        f.write('_atom_site_occupancy\n')

        for i in range(len(atoms['position'])):
            f.write("{:14s}{:4s}{:12.6f}{:12.6f}{:12.6f}   1.0000\n".format(
                atoms['label'][i],atoms['symbol'][i],atoms['position'][i][0],
                atoms['position'][i][1],atoms['position'][i][2]))

        if (not np.isnan(atoms['Uij']).any()):
            f.write("loop_\n")
            f.write("_atom_site_aniso_label\n")
            f.write("_atom_site_aniso_U_11\n")
            f.write("_atom_site_aniso_U_22\n")
            f.write("_atom_site_aniso_U_33\n")
            f.write("_atom_site_aniso_U_23\n")
            f.write("_atom_site_aniso_U_13\n")
            f.write("_atom_site_aniso_U_12\n")
            for i in range(len(atoms['label'])):
                f.write("{:14s} {:9.5f}{:9.5f}{:9.5f}{:9.5f}{:9.5f}{:9.5f}\n".format(
                    atoms['label'][i],
                    atoms['Uij'][i][0][0],
                    atoms['Uij'][i][1][1],
                    atoms['Uij'][i][2][2],
                    atoms['Uij'][i][1][2],
                    atoms['Uij'][i][0][2],
                    atoms['Uij'][i][0][1]))

        if verbose:
            print('data_cif\n')
            print('{:32s}{:7.4f}(0)'.format('_cell_length_a',lat['L_a']))
            print('{:32s}{:7.4f}(0)'.format('_cell_length_b',lat['L_b']))
            print('{:32s}{:7.4f}(0)'.format('_cell_length_c',lat['L_c']))
            print('{:32s}{:7.4f}(0)'.format('_cell_angle_alpha',degrees(lat['alpha'])))
            print('{:32s}{:7.4f}(0)'.format('_cell_angle_beta',degrees(lat['beta'])))
            print('{:32s}{:7.4f}(0)\n'.format('_cell_angle_gamma',degrees(lat['gamma'])))
            print("{:35s} 'P1'".format('_symmetry_space_group_name_H-M'))
            print('{:35s} 1'.format('_symmetry_Int_Tables_number'))
            print('{:35s} triclinic\n'.format('_symmetry_cell_setting'))
            print('loop_')
            print('_atom_site_label')
            print('_atom_site_type_symbol')
            print('_atom_site_fract_x')
            print('_atom_site_fract_y')
            print('_atom_site_fract_z')
            print('_atom_site_occupancy')

            for i in range(len(atoms['position'])):
                print("{:14s}{:4s}{:12.6f}{:12.6f}{:12.6f}   1.0000".format(
                    atoms['label'][i],atoms['symbol'][i],atoms['position'][i][0],
                    atoms['position'][i][1],atoms['position'][i][2]))

            if (not np.isnan(atoms['Uij']).any()):
                print(
                    "loop_\n",
                    "_atom_site_aniso_label\n",
                    "_atom_site_aniso_U_11\n",
                    "_atom_site_aniso_U_22\n",
                    "_atom_site_aniso_U_33\n",
                    "_atom_site_aniso_U_23\n",
                    "_atom_site_aniso_U_13\n",
                    "_atom_site_aniso_U_12")
                for i in range(len(atoms['label'])):
                    print("{:14s} {:9.5f}{:9.5f}{:9.5f}{:9.5f}{:9.5f}{:9.5f}".format(
                        atoms['label'][i],
                        atoms['Uij'][i][0][0],
                        atoms['Uij'][i][1][1],
                        atoms['Uij'][i][2][2],
                        atoms['Uij'][i][1][2],
                        atoms['Uij'][i][0][2],
                        atoms['Uij'][i][0][1]))

def split_cif(filename,dir='./',name_delim=None): 
    '''
    Routine for reading .CIF files that contain multiple 
    Places cif file in new style dictionary and an Xtal class, then writes to directory with data name
    Presently requires plain text descriptions of symetry ops, and does not interpret space groups
    ###############
    Parameters:
        filename: complete location of cif file
        dir:  target directory
        name_delim: splits data_ string with deliminator and rejoins with underscore '_'
    ###############
    '''
    #===========================================================================
    # # FOR CHENGXI'S DATA ONLY #
    # import pandas as pd
    # df = pd.DataFrame()
    # df = pd.DataFrame(columns=['Name', 'Density', 'Energy'])
    # basename = os.path.splitext(os.path.basename(filename))[0]
    # # FOR CHENGXI'S DATA ONLY #
    #===========================================================================
    # Open the .CIF file
    with open(filename,'r') as f:
        line_num = 0
        cif_num = -1
        same_cif = True
        first_cif = True
        reading_sym_ops = False
        reading_atom_sites = False
        reading_atom_aniso = False
        data={}
        cell_data={}
        # Read lines one by one
        for line in f:
            line_num += 1
            if not same_cif:
                reading_sym_ops = False
                reading_atom_sites = False
                reading_atom_aniso = False
                same_cif = True
                data={}
                cell_data={}         
                
            # Split into columns
            cols = line.split()
            if (len(cols) == 0): continue
            # Pull data title
            if cols[0][:5]=='data_':
                cif_num+=1
                data_name_tmp = cols[0][5:]
                if name_delim:
                    data_name_tmp = '_'.join(data_name_tmp.split(name_delim))
                #===============================================================
                # # FOR CHENGXI'S DATA ONLY #
                # cols = cols[0].split('_')
                # energy = float(cols[2])
                # density = float(cols[3])
                # df.loc[cif_num] = ["{}_{}".format(basename,cif_num),density,energy] 
                # # FOR CHENGXI'S DATA ONLY #
                #===============================================================
                if first_cif: 
                    first_cif=False
                    data_name=data_name_tmp
                    continue
            
            # Output cif at end of record
            if cols[0][:4]=='#END' or cols[0][:5]=='data_':
                data_name_tmp = cols[0][5:]
                if name_delim:
                    data_name_tmp = '_'.join(data_name_tmp.split(name_delim))
                same_cif=False
                #Creating a cleaner dictionary to return with a complete atom list
                cell_data['L_a'] =  float(data['_cell_length_a'])
                cell_data['L_b'] = float(data['_cell_length_b'])
                cell_data['L_c'] = float(data['_cell_length_c'])
                cell_data['alpha'] = radians(float(data['_cell_angle_alpha']))
                cell_data['beta'] = radians(float(data['_cell_angle_beta']))
                cell_data['gamma'] = radians(float(data['_cell_angle_gamma']))
                if ('_cell_volume' in data): cell_data['volume'] = float(data['_cell_volume'])
            
                #Creating atom list
                syms = data['_atom_site_type_symbol']
                labels = data['_atom_site_label']
                coords = np.transpose(np.array([data['_atom_site_fract_x'],data['_atom_site_fract_y'],
                                   data['_atom_site_fract_z']]))
                
                if ('_atom_site_aniso_label' in data):
                    Uij = np.transpose(np.array([
                            [data['_atom_site_aniso_U_11'],data['_atom_site_aniso_U_12'],data['_atom_site_aniso_U_13']],
                            [data['_atom_site_aniso_U_12'],data['_atom_site_aniso_U_22'],data['_atom_site_aniso_U_23']],
                            [data['_atom_site_aniso_U_13'],data['_atom_site_aniso_U_23'],data['_atom_site_aniso_U_33']]
                                    ]))
                else:
                    Uij=None
            
                # Fixing with float division and applying symmetry operations
                if ('_symmetry_equiv_pos_as_xyz' in data):
                    ops = data['_symmetry_equiv_pos_as_xyz']
                    for i in range(len(ops)):
                        ops[i] = ops[i].replace("/","./")
            
                    # Two atoms are on top of each other if they are less than "eps" away.
                    eps = 0.01  # in Angstrom
                    imax = len(syms)
                    new_coords=[]
                    new_syms=[]
                    new_labels=[]
                    new_U=[]
                    for i in range(imax):
                        for op in ops:
                            x,y,z = coords[i,:]
                            # Text evaluation of symmetry operation
                            xn = numexpr.evaluate(op.split(',')[0]).item()
                            yn = numexpr.evaluate(op.split(',')[1]).item()
                            zn = numexpr.evaluate(op.split(',')[2]).item()
                            # Forcing into the unit cell
                            xn = (xn + 10.0) % 1.0
                            yn = (yn + 10.0) % 1.0
                            zn = (zn + 10.0) % 1.0
                            # Adding only unique new atoms
                            new_atom = True
                            for c in range(imax):
                                if (abs(coords[c,0]-xn)<eps and abs(coords[c,1]-yn)<eps and abs(coords[c,2]-zn)<eps):
                                    new_atom = False
                                    break
                            for c in range(len(new_coords)):
                                if (abs(new_coords[c][0]-xn)<eps and abs(new_coords[c][1]-yn)<eps and abs(new_coords[c][2]-zn)<eps):
                                    new_atom = False
                                    break
                            if new_atom:
                                new_coords.append(np.array([xn,yn,zn]))
                                new_syms.append(syms[i])
                                new_labels.append(labels[i])
                                if Uij: new_U.append(Uij[i,:])
                
                #Sorting the atom list alphabetically
                if (new_coords):
                    coords=np.vstack((coords, new_coords))
                    syms.extend(new_syms)
                    labels.extend(new_labels)
                    if Uij: Uij=np.vstack((Uij,new_U))
                xtal = Xtal(positions=coords,lat=cell_data,symbols=syms,labels=labels,Uij=Uij)
                write_file = os.path.join(dir,"{}.cif".format(data_name))
                print(write_file)
                write_cif(write_file,xtal) 
                data_name=data_name_tmp
                continue
                    
            # ID the keywords
            # Keywords associated with cell parameters
            if (cols[0] == '_cell_length_a'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_a'] = float(cols[1])
            elif (cols[0] == '_cell_length_b'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_b'] = float(cols[1])
            elif (cols[0] == '_cell_length_c'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_c'] = float(cols[1])
            elif (cols[0] == '_cell_angle_alpha'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_alpha'] = float(cols[1])
            elif (cols[0] == '_cell_angle_beta'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_beta'] = float(cols[1])
            elif (cols[0] == '_cell_angle_gamma'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_gamma'] = float(cols[1])
            elif (cols[0] == '_cell_volume'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_volume'] = float(cols[1])

            # Keywords associated with symmetry operations
            elif (cols[0] == '_symmetry_equiv_pos_as_xyz'):
                reading_sym_ops = True
                data['_symmetry_equiv_pos_as_xyz'] = []
            elif (reading_sym_ops):
                # Add the operation if the string is between single quotes.
                # Otherwise it's a sign we are done with the list.
                # Some files don't use the quotes, so wait for next underscored command or loop_
                if (cols[0][0] == '\''  and  cols[0][-1] == '\''):
                    data['_symmetry_equiv_pos_as_xyz'].append(cols[0][1:-1])
                elif cols[0][0] == '_' or cols[0][:5]=='loop_':
                    reading_sym_ops = False
                else:
                    data['_symmetry_equiv_pos_as_xyz'].append(cols[0][:])

            # Keywords associated with atom sites
            elif (cols[0] == '_atom_site_label'):
                index_label = line_num
                data['_atom_site_label'] = []
                data['_atom_site_fract_x'] = []
                data['_atom_site_fract_y'] = []
                data['_atom_site_fract_z'] = []
                reading_atom_sites = True

            # Keep track of where the other labels are (order is important).
            elif (cols[0] == '_atom_site_fract_x'):
                index_x = line_num - index_label
            elif (cols[0] == '_atom_site_fract_y'):
                index_y = line_num - index_label
            elif (cols[0] == '_atom_site_fract_z'):
                index_z = line_num - index_label
            elif (cols[0] == '_atom_site_type_symbol'):
                index_sym = line_num - index_label
                data['_atom_site_type_symbol'] = []

            # If we are currently reading the atom sites...
            elif (reading_atom_sites):
                # Read the actual data if we have 4 columns or more of data.
                if (len(cols) >= 4):
                    data['_atom_site_label'].append(cols[0])
                    data['_atom_site_fract_x'].append(float(cols[index_x].split('(',1)[0]))
                    data['_atom_site_fract_y'].append(float(cols[index_y].split('(',1)[0]))
                    data['_atom_site_fract_z'].append(float(cols[index_z].split('(',1)[0]))
                    if ('_atom_site_type_symbol' in data): data['_atom_site_type_symbol'].append(cols[index_sym])

                # Stop reading atom sites if we found a line with fewer
                # columns, and which does not start with '_atom_site_'.
                elif (len(cols[0]) < 11  or  cols[0][:11] != '_atom_site_'):
                    reading_atom_sites = False
            
            #Keywords associated with aniso. Presently only works for isoptropic case with symmetry. Anisotropic without
            elif(cols[0] == '_atom_site_aniso_label'):
                index_label = line_num
                data['_atom_site_aniso_label'] = []
                data['_atom_site_aniso_U_11'] = []
                data['_atom_site_aniso_U_22'] = []
                data['_atom_site_aniso_U_33'] = []
                data['_atom_site_aniso_U_23'] = []
                data['_atom_site_aniso_U_13'] = []
                data['_atom_site_aniso_U_12'] = []
                reading_atom_aniso = True

            # Keep track of where the other labels are (order is important).
            elif(cols[0] == '_atom_site_aniso_U_11'):
                index_11 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_22'):
                index_22 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_33'):
                index_33 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_23'):
                index_23 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_13'):
                index_13 = line_num-index_label
            elif(cols[0] == '_atom_site_aniso_U_12'):
                index_12 = line_num-index_label
            elif(reading_atom_aniso):
                # Read the actual data if we have 4 columns or more of data.
                if (len(cols) >= 4):
                    data['_atom_site_aniso_label'].append(cols[0])
                    data['_atom_site_aniso_U_11'].append(float(cols[index_11]))
                    data['_atom_site_aniso_U_22'].append(float(cols[index_22]))
                    data['_atom_site_aniso_U_33'].append(float(cols[index_33]))
                    data['_atom_site_aniso_U_23'].append(float(cols[index_23]))
                    data['_atom_site_aniso_U_13'].append(float(cols[index_13]))
                    data['_atom_site_aniso_U_12'].append(float(cols[index_12]))

                # Stop reading atom sites if we found a line with fewer
                # columns, and which does not start with '_atom_site_'.
                elif (len(cols[0]) < 11  or  cols[0][:11] != '_atom_site_'):
                    reading_atom_aniso = False
    #===========================================================================
    # # FOR CHENGXI'S DATA ONLY #
    # df.to_pickle(os.path.join(dir,"{}_CSP_data".format(basename)))
    #===========================================================================
    