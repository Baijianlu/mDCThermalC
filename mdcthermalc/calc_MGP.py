# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 11:30:45 2019

@author: Tao.Fan
This script for calculating mode Gruneisen by averaging points around Gamma points
"""
import numpy as np
from mdcthermalc.extract_GrunP import extract_GrunP
from mdcthermalc.extract_GV import extract_GV

def calc_MGP(filepath,weight):                        #Gamma:the position of Gamma, weight:multiplicity of Gamma points
    (GroupVec,Freq,Gamma) = extract_GV(filepath)
#    Gamma = np.array([0.1363090,0.5483771,0.9786495])              #Gamma point should be decided by external tools
    Gamma_index = np.zeros((len(Gamma),2))          #Gamma point maybe the first one, thus the index could be zero, must be careful
    Bandnum = GroupVec.shape[1] - 1                 ###All these parameter can be transfered into function in the future 
    Minindex = np.int(0)
    Maxindex = np.int(GroupVec.shape[0] - 1)    
    Path, Gruneisen = extract_GrunP(filepath,Bandnum,GroupVec.shape[0])
#    (no1,no2,weight) = get_highsympath("Mg2Si_mp-1367.cif")
    #搜索groupvec得到Gamma点的位置，并向gamma点左右两边各取5个点
    Pathnum = 0
    for i in np.arange(len(Gamma)):
        Gamma_index[i] = np.array([x_index for x_index,x_value in enumerate(GroupVec[:,0]) if x_value==Gamma[i]])
        if Gamma_index[i,0] == Gamma_index[i,1]:Pathnum = Pathnum + 1
        else: Pathnum = Pathnum + 2
        
    modebranch_grun = np.zeros((Bandnum,Pathnum,51))        # value for each path in each branch
    branch_grun = np.zeros(Bandnum)             #average value for different branch
    
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

    
#    [np.savetxt('modegruneisen' + np.str(i) + '.out',modebranch_grun[i],fmt='%.8f',delimiter=',') for i in np.arange(Bandnum)]
    return branch_grun

if __name__ == '__main__':
   branch_grun = calc_MGP()
