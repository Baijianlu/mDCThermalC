# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 15:16:02 2019

@author: Tao.Fan 
this script for calculating the gruneisen,velocity,DebyeT used for kappa calculation,
they are all four dimension
"""
import numpy as np
import scipy.constants
from mdcthermalc.calc_MGP import calc_MGP
from mdcthermalc.calc_MGV import calc_MGV
from mdcthermalc.get_highsymweight import get_highsymweight

def Get_GVD(filepath):
    planck = scipy.constants.h
    Boltzm = scipy.constants.Boltzmann
    gruneisen = np.zeros(4)
    velocity = np.zeros(4)
    DebyeT = np.zeros(4)
    freq = np.zeros(4)
    
    (no1,no2,weight) = get_highsymweight(filepath + "/POSCAR")
    (branchvel,branchfreq,branchDebyeT) = calc_MGV(filepath,weight)
    branchgrun = calc_MGP(filepath,weight)
    
    gruneisen[0:3] = branchgrun[0:3]
    velocity[0:3] = branchvel[0:3]
    DebyeT[0:3] = branchDebyeT[0:3]
    freq[0:3] = branchfreq[0:3]
    weightsum = np.sum(branchfreq[3:])
    
    for i in np.arange(3,len(branchfreq)):
        gruneisen[3] = gruneisen[3] + branchfreq[i] * branchgrun[i]
        velocity[3] = velocity[3] + branchfreq[i] * branchvel[i]
        DebyeT[3] = DebyeT[3] + branchfreq[i] * branchDebyeT[i]
    
    gruneisen[3] = gruneisen[3]/weightsum
    velocity[3] = velocity[3]/weightsum
    DebyeT[3] = DebyeT[3]/weightsum
    freq[3] = DebyeT[3] * Boltzm/(1e12 * planck)
    return gruneisen, velocity, DebyeT, freq

if __name__ == '__main__':
   (gruneisen,velocity,DebyeT,freq) = Get_GVD()
