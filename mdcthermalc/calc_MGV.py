# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 12:43:43 2019

@author: Tao.Fan
This script for calculating mode group velocity by averaging points around Gamma points
"""
import numpy as np
import scipy.constants
from mdcthermalc.extract_GV import extract_GV

def calc_MGV(filepath,weight):
    planck = scipy.constants.h
    Boltzm = scipy.constants.Boltzmann
    
    (GroupVec,Frequency,Gamma) = extract_GV(filepath)
#    Gamma = np.array([0.1363090,0.5483771,0.9786495])              #Gamma point should be decided by external tools   Mg2Si:0.1363090,0.5483771   LiMg[0.14588740,0.75077500]   AL3Li[0.12431940,0.63978040]
    Gamma_index = np.zeros((len(Gamma),2))          #Gamma point maybe the first one, thus the index could be zero, must be careful
    Bandnum = GroupVec.shape[1] - 1
    Minindex = np.int(0)
    Maxindex = np.int(GroupVec.shape[0] - 1)
#    weight = np.array([8,12,48,24,24])                  #this should be decided by function   Mg2Si:[8,6,24,12]   LiMg:[6,12,8] Al3Li:[6,12,8]
#    (no1,no2,weight) = get_highsympath("Mg2Si_mp-1367.cif")
    #搜索groupvec得到Gamma点的位置，并向gamma点左右两边各取5个点
    Pathnum = 0
    for i in np.arange(len(Gamma)):
        Gamma_index[i] = np.array([x_index for x_index,x_value in enumerate(GroupVec[:,0]) if x_value==Gamma[i]])
    ##此处需要有个语句来判断有多少条路径出现，从而确定modebranch_vel的第二维的大小    
        if Gamma_index[i,0] == Gamma_index[i,1]:Pathnum = Pathnum + 1
        else: Pathnum = Pathnum + 2
    
    ####
    modebranch_vel = np.zeros((Bandnum,Pathnum,5))   #the first dimension size equal natom*3
    branch_vel = np.zeros(Bandnum)
    
    for branch_idx in np.arange(Bandnum):
        for j in np.arange(len(Gamma)):                 #此处j的范围是非对称路径的条数
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
    
#    [np.savetxt('modevelocity' + np.str(i) + '.out',modebranch_vel[i],fmt='%.8f',delimiter=',') for i in np.arange(Bandnum)]
    
    ###the following is for calculating average frequency for different branch
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
    
    ###the following is for calculating debye temperature for different branch
    branch_DebyeT = np.zeros(Bandnum)
    for branch_idx in np.arange(Bandnum):
        branch_DebyeT[branch_idx] = planck * np.max(Frequency[:,branch_idx+1]) * 1e12/Boltzm
    
    return branch_vel,branch_freq,branch_DebyeT


if __name__ == '__main__':
   (branchv,branchf,branchD) = calc_MGV()
