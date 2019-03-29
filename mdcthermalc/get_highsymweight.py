# -*- coding: utf-8 -*-
"""
Copyright (C) 2019 Tao Fan 
All rights reserved

This script used to calculate the weight of each path of high symmetry paths. Such a weight is equal to symmetrical equivalent point
number of each point in that path.

"""
import numpy as np
import pymatgen as pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

def Coordcharacter(coord):                                         #coord is a np.array format
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

def get_highsymweight(filename):
    Mg2Si = pmg.Structure.from_file(filename)       
    finder = SpacegroupAnalyzer(Mg2Si)
    symbol = finder.get_space_group_symbol()
    HKpath = HighSymmKpath(Mg2Si)
    Keys = list()
    Coords = list()
    
    for key in HKpath.kpath['kpoints']:
        Keys.append(key)
        Coords.append(HKpath.kpath['kpoints'][key])
        
    count = 0
    Keylist = list()
    Coordslist = list()
    for i in np.arange(len(Keys) - 1):
        if (count-1)%3 == 0:                                        #count-1 can be intergely divided by 3
            Keylist.append(Keys[0])
            Coordslist.append(Coords[0])
            count+=1
            
        Keylist.append(Keys[i+1])
        Coordslist.append(Coords[i+1])
        count+=1
        
    if (count-1)%3 == 0:
        Keylist.append(Keys[0])
        Coordslist.append(Coords[0])
        
    Kweight = np.zeros(len(Keys) - 1)
    kmesh = finder.get_ir_reciprocal_mesh(mesh=(50,50,50))
        
    for i in np.arange(len(Keys) - 1):
        (zerocount,nonzeroratio) = Coordcharacter(Coords[i+1])
        
        for j in np.arange(len(kmesh)):
            (mzerocount,mnonzeroratio) = Coordcharacter(kmesh[j][0])
            if len(mnonzeroratio) == len(nonzeroratio):
                remainlogic = np.abs(nonzeroratio - mnonzeroratio) < 0.01                  # 0.01 is a value can get enough accurate results
                if zerocount == mzerocount and remainlogic.all():
                    if kmesh[j][1] > Kweight[i]:
                        Kweight[i] = kmesh[j][1]
            else: pass
    
    return Keylist, Coordslist, Kweight


if __name__ == '__main__':
   (Keylist,Coordslist,Kweight) = get_highsymweight("POSCAR")

                
