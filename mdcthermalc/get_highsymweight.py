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

def pbc_diff(fcoords1, fcoords2):
    fdist = np.subtract(fcoords1, fcoords2)
    return fdist - np.round(fdist)
    
def get_sym_eq_kpoints(struct, kpoint, cartesian=False, tol=1e-2):
    if not struct:
        return None
    
    sg = SpacegroupAnalyzer(struct)
    symmops = sg.get_point_group_operations(cartesian=cartesian)
    points = np.dot(kpoint, [m.rotation_matrix for m in symmops])
    rm_list = []
    # identify and remove duplicates from the list of equivalent k-points:
    for i in range(len(points) - 1):
        for j in range(i + 1, len(points)):
            if np.allclose(pbc_diff(points[i], points[j]), [0, 0, 0], tol):
                rm_list.append(i)
                break
    return np.delete(points, rm_list, axis=0)

def get_highsymweight(filename):
    struct = pmg.Structure.from_file(filename)       #here should be changed
    HKpath = HighSymmKpath(struct)
    Keys = list()
    Coords = list()
        
    for key in HKpath.kpath['kpoints']:
        Keys.append(key)
        Coords.append(HKpath.kpath['kpoints'][key])
    
    Kweight = list()
        
    for i in np.arange(len(Keys)):
        if Keys[i] != '\Gamma':
           Kweight.append(len(get_sym_eq_kpoints(struct, Coords[i]*0.5))) 
        
    return Keys, Coords, Kweight


if __name__ == '__main__':
   (Keylist,Coordslist,Kweight) = get_highsymweight("POSCAR")

                
