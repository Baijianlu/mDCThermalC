# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 15:05:43 2019

@author: Tao.Fan
This file for extracting group velocity information in band.yaml

"""
import os
import numpy as np

def extract_GV(filepath):
    fp1 = open(filepath + '/band.yaml','r')

    keystr1 = "q-position"
    keystr2 = "distance"
    keystr3 = "frequency"
    keystr4 = "group_velocity"
    keystr5 = "nqpoint:"
    keystr6 = "natom:"
    npoints = 0
    nbands = 0
    countpoints = -1
    countbands = 0
    Gammaflag = 0
    Gamma = list()

    for eachline in fp1:
        eachline = eachline.strip()
        temp = eachline.split()
        if len(temp) > 0:
            if keystr5 == temp[0]:        
                npoints = int(temp[-1])
            elif keystr6 in eachline:
                nbands = int(temp[-1]) * 3
                GroupVec = np.zeros((npoints, nbands+1))
                Frequency = np.zeros((npoints, nbands+1))
            elif keystr1 in eachline:
                countpoints = countpoints + 1
                countbands = 0
                postemp = np.array([np.float(temp[i][:-1]) for i in np.arange(3,6)])
#                print('%f %f %f' % (postemp[0],postemp[1],postemp[2]))
                if postemp[0] == 0.0 and postemp[1] == 0.0 and postemp[2] == 0.0:
                    Gammaflag = 1
            elif keystr2 in eachline:
                #write distance value to the first column of np.array
                GroupVec[countpoints,countbands] = np.float(temp[-1])
                Frequency[countpoints,countbands] = np.float(temp[-1])
                countbands = countbands + 1
                if Gammaflag == 1:
                    Gammaflag = 0
                    if np.float(temp[-1]) not in Gamma:                
                        Gamma.append(np.float(temp[-1]))
            elif keystr3 in eachline:
                Frequency[countpoints,countbands] = np.float(temp[-1])
            elif keystr4 in eachline:
                #write velocity value to the rest colume of each row of np .array
                vectemp = np.array([np.float(temp[i][:-1]) for i in np.arange(2,5)])
                vectemp2 = vectemp**2
                GroupVec[countpoints,countbands] = np.sqrt(vectemp2.sum())
                countbands = countbands + 1
            else:
                continue
        else:
            continue

    fp1.close()
    Gamma = np.array(Gamma)

    return GroupVec,Frequency,Gamma

if __name__ == '__main__':
   (GroupV,Freqc,Gamma) = extract_GV()




