# -*- coding: utf-8 -*-
"""
Copyright (C) 2019 Tao Fan 
All rights reserved

This script is used for extracting Gruneisen parameter in gruneisen.yaml 

"""
import numpy as np

def extract_GrunP(filepath, nbands=9, npoints=255):
    
    fp1 = open(filepath + '/gruneisen.yaml','r')

    keystr1 = "q-position"
    keystr2 = "distance"
    keystr3 = "gruneisen"
    keystr4 = "frequency"
    datatype = np.dtype([('freq',np.float),('grun',np.float)])
    Pathpot = np.zeros(npoints)
    GruneisenPara = np.zeros((npoints, nbands),dtype=datatype)
    countpoints = -1
    countbands = 0
    
    
    for eachline in fp1:
        eachline = eachline.strip()
        temp = eachline.split()
        if len(temp) > 0:
            if keystr1 in eachline:
                countpoints = countpoints + 1
                countbands = 0
            elif keystr2 in eachline:
                #write distance value to the first column of np.array
                Pathpot[countpoints] = np.float(temp[-1])            
            elif keystr3 in eachline:
                #write gruneisen value to the rest colume of each row of np.array
                gruntemp = np.float(temp[-1])            
            elif keystr4 in eachline:
                freqtemp = np.float(temp[-1])
                GruneisenPara[countpoints,countbands] = np.array((freqtemp,gruntemp), dtype=datatype)
                countbands = countbands + 1
            else:
                continue
        else:
            continue
    
    
    fp1.close()
    
    for j in range(npoints):
        GruneisenPara[j,:] = np.sort(GruneisenPara[j,:], order='freq')
        
    return Pathpot, GruneisenPara    

if __name__ == '__main__':
   Path, Grunsisen = extract_GrunP()
