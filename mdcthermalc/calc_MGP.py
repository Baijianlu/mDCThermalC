# -*- coding: utf-8 -*-
"""
Copyright (C) 2019 Tao Fan 
All rights reserved

This script is used for calculating Gruneisen parameter of each branch by weight-averaged Gruneisen parameter of each path in the branch.

"""
import numpy as np
from mdcthermalc.extract_GrunP import extract_GrunP
from mdcthermalc.extract_GV import extract_GV

def calc_MGP(filepath,weight):                                     #Gamma:the position of Gamma, weight:multiplicity of Gamma points
    (GroupVec,Freq,Gamma) = extract_GV(filepath)
    Gamma_index = np.zeros((len(Gamma),2))          
    Bandnum = GroupVec.shape[1] - 1                 
    Minindex = np.int(0)
    Maxindex = np.int(GroupVec.shape[0] - 1)    
    Path, Gruneisen = extract_GrunP(filepath,Bandnum,GroupVec.shape[0])
    Pathnum = 0
    
    for i in np.arange(len(Gamma)):
        Gamma_index[i] = np.array([x_index for x_index,x_value in enumerate(GroupVec[:,0]) if x_value==Gamma[i]])
        if Gamma_index[i,0] == Gamma_index[i,1]:Pathnum = Pathnum + 1
        else: Pathnum = Pathnum + 2
        
    modebranch_grun = np.zeros((Bandnum,Pathnum,51))              # value for each path in each branch
    branch_grun = np.zeros(Bandnum)                               #average value for different branch
    
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(len(Gamma)): 
            for k in np.arange(2):
                if k == 0:
                    if Gamma_index[j,k] > Minindex:
                        modebranch_grun[branch_idx,j*2 + k] = [Gruneisen[np.int(index),branch_idx]['grun'] for index in np.arange(Gamma_index[j,k]-50,Gamma_index[j,k]+1)]
                    else:
                        modebranch_grun[branch_idx,j*2 + k] = [Gruneisen[np.int(index),branch_idx]['grun'] for index in np.arange(Gamma_index[j,k],Gamma_index[j,k]+51)]
                        break
                if k == 1:
                    if Gamma_index[j,k] < Maxindex:
                        modebranch_grun[branch_idx,j*2 + k] = [Gruneisen[np.int(index),branch_idx]['grun'] for index in np.arange(Gamma_index[j,k],Gamma_index[j,k]+51)]
                    else:
                        break
                    
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(len(weight)):
            branch_grun[branch_idx] = branch_grun[branch_idx] + weight[j] * np.power(np.average(modebranch_grun[branch_idx,j,:]),2)
        branch_grun[branch_idx] = np.sqrt(branch_grun[branch_idx] / np.sum(weight))    

    
    return branch_grun

if __name__ == '__main__':
   branch_grun = calc_MGP()
