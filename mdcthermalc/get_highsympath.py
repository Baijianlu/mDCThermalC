# -*- coding: utf-8 -*-
"""
Copyright (C) 2019 Tao Fan 
All rights reserved

This script used to obtain the exact high symmetry path that should be used for phonon spectrum calculation. The high symmetry path are 
aranged such that each high symmetry point beside Gamma point have one and only have one connection with Gamma point.

"""
import os
import numpy as np
from numpy.linalg import inv
import pymatgen as pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

def Coordcharacter(coord):                                      #coord is a np.array format
    zerocount = np.int(0)
    nonzeroratio = list()
    
    nonPos = coord.nonzero()[0]
    zerocount = 3 - len(nonPos)
    
    if zerocount == 0:
        nonzeroratio.append(coord[0]/coord[1])
        nonzeroratio.append(coord[1]/coord[0])
        nonzeroratio.append(coord[0]/coord[2])
        nonzeroratio.append(coord[2]/coord[0])
        nonzeroratio.append(coord[1]/coord[2])
        nonzeroratio.append(coord[2]/coord[1])
    elif zerocount == 1:
        nonzeroratio.append(coord[nonPos[0]]/coord[nonPos[1]])
        nonzeroratio.append(coord[nonPos[1]]/coord[nonPos[0]])
    elif zerocount == 2:
        nonzeroratio.append(coord[nonPos[0]]/coord[nonPos[0]])
    else:
        pass
    
    nonzeroratio.sort()
    
    return zerocount, np.array(nonzeroratio)

def get_highsympath(filename):
    struct = pmg.Structure.from_file(filename)      
    finder = SpacegroupAnalyzer(struct)
    HKpath = HighSymmKpath(struct)
    Keys = list()
    Coords = list()
    
    for key in HKpath.kpath['kpoints']:
        Keys.append(key)
        Coords.append(HKpath.kpath['kpoints'][key])
        
    count = 0
    Keylist = list()
    Coordslist = list()
    for i in np.arange(len(Keys) - 1):
        if (count-1)%3 == 0:                                   #count-1 can be intergely divided by 3
            Keylist.append(Keys[0])
            Coordslist.append(Coords[0])
            count+=1
            
        Keylist.append(Keys[i+1])
        Coordslist.append(Coords[i+1])
        count+=1
        
    if (count-1)%3 == 0:
        Keylist.append(Keys[0])
        Coordslist.append(Coords[0])
        
    prims = finder.find_primitive()
    
    print('Please set \"BAND\" parameter of phonopy as this:%s' % os.linesep)
    for coord in Coordslist:
        print('%.4f %.4f %.4f  ' % (coord[0], coord[1], coord[2]), end='')
    print('%s' % os.linesep)
    
    if struct != prims:
        S_T = np.transpose(struct.lattice.matrix)
        P_T = np.transpose(prims.lattice.matrix)
        transmat = inv(S_T) @ P_T
        print('We notice your structure could have a primitive cell. Please set \"PRIMITIVE_AXIS\" parameter of phonopy as this:%s' % os.linesep)
        for coord in transmat:
            print('%.8f %.8f %.8f  ' % (coord[0], coord[1], coord[2]), end='')
        print('%s' % os.linesep)
        
        
    return Keylist, Coordslist


if __name__ == '__main__':
   (Keylist,Coordslist) = get_highsympath("POSCAR")
