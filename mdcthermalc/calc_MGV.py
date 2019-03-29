# -*- coding: utf-8 -*-
"""
Copyright (C) 2019 Tao Fan 
All rights reserved

This script is used for calculating group velocity, Debye temperature of each branch. The group velocity of branch is weighted average value of each path 
in this branch. Debye temperature of each branch is calculated by using the highest frequency of the branch. 

"""
import numpy as np
import scipy.constants
from mdcthermalc.extract_GV import extract_GV

def calc_MGV(filepath,weight):
    planck = scipy.constants.h
    Boltzm = scipy.constants.Boltzmann
    
    (GroupVec,Frequency,Gamma) = extract_GV(filepath)
    Gamma_index = np.zeros((len(Gamma),2))          
    Bandnum = GroupVec.shape[1] - 1
    Minindex = np.int(0)
    Maxindex = np.int(GroupVec.shape[0] - 1)

    #search groupvec to get Gamma positions
    Pathnum = 0
    for i in np.arange(len(Gamma)):
        Gamma_index[i] = np.array([x_index for x_index,x_value in enumerate(GroupVec[:,0]) if x_value==Gamma[i]])   
        if Gamma_index[i,0] == Gamma_index[i,1]:Pathnum = Pathnum + 1
        else: Pathnum = Pathnum + 2
    
    #the following is for calculating average group velocity of different branch
    modebranch_vel = np.zeros((Bandnum,Pathnum,5))   #the first dimension size equal natom*3
    branch_vel = np.zeros(Bandnum)
    
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(len(Gamma)):                 #
            for k in np.arange(2):
                if k == 0: 
                    if Gamma_index[j,k] > Minindex:
                        modebranch_vel[branch_idx,j*2 + k] = [GroupVec[np.int(index),branch_idx+1] for index in np.arange(Gamma_index[j,k]-5,Gamma_index[j,k])]
                    else:                            #actually, this sentence should never be executed and G point should never be the first point
                        modebranch_vel[branch_idx,j*2 + k] = [GroupVec[np.int(index),branch_idx+1] for index in np.arange(Gamma_index[j,k]+1,Gamma_index[j,k]+6)]
                        break
                if k == 1: 
                    if Gamma_index[j,k] < Maxindex:
                        modebranch_vel[branch_idx,j*2 + k] = [GroupVec[np.int(index),branch_idx+1] for index in np.arange(Gamma_index[j,k]+1,Gamma_index[j,k]+6)]
                    else:
                        break
                
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(Pathnum):
            branch_vel[branch_idx] = branch_vel[branch_idx] + weight[j] * np.average(modebranch_vel[branch_idx,j,:])
        branch_vel[branch_idx] = branch_vel[branch_idx] / np.sum(weight)
    
    #the following is for calculating average frequency of different branch
    modebranch_freq = np.zeros((Bandnum,Pathnum,51))
    branch_freq = np.zeros(Bandnum)
    
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(len(Gamma)): 
            for k in np.arange(2):
                if k == 0:
                    if Gamma_index[j,k] > Minindex:
                        modebranch_freq[branch_idx,j*2 + k] = [Frequency[np.int(index),branch_idx+1] for index in np.arange(Gamma_index[j,k]-50,Gamma_index[j,k]+1)]
                    else:
                        modebranch_freq[branch_idx,j*2 + k] = [Frequency[np.int(index),branch_idx+1] for index in np.arange(Gamma_index[j,k],Gamma_index[j,k]+51)]
                        break
                if k == 1:
                    if Gamma_index[j,k] < Maxindex:
                        modebranch_freq[branch_idx,j*2 + k] = [Frequency[np.int(index),branch_idx+1] for index in np.arange(Gamma_index[j,k],Gamma_index[j,k]+51)]
                    else:
                        break
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(Pathnum):
            branch_freq[branch_idx] = branch_freq[branch_idx] + weight[j] * np.average(modebranch_freq[branch_idx,j,:])
        branch_freq[branch_idx] = branch_freq[branch_idx] / np.sum(weight)    
    
    #the following is for calculating debye temperature for different branch
    branch_DebyeT = np.zeros(Bandnum)
    for branch_idx in np.arange(Bandnum):
        branch_DebyeT[branch_idx] = planck * np.max(Frequency[:,branch_idx+1]) * 1e12/Boltzm
    
    return branch_vel,branch_freq,branch_DebyeT


if __name__ == '__main__':
   (branchv,branchf,branchD) = calc_MGV()
