# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 15:11:56 2019
@author: Tao.Fan
This script for calculating mass-fluctuation phonon-scattering parameter Gamma
"""
import numpy as np
import pymatgen as pmg

MassFluct = {'H':1.1460e-4, 'He':8.3232e-8, 'Li':14.58e-4, 'Be':0.0, 'B':13.54e-4, 'C':7.38695e-05, 'N':1.8577e-05, 'O':3.3590e-05, 'F':0.0, 'Ne':8.2792e-4, 
             'Na':0.0, 'Mg':7.3989e-4, 'Al':0.0, 'Si':2.01222e-4, 'P':0.0, 'S':1.6808e-4, 'Cl':5.8237e-4, 'Ar':3.50987e-05,
             'K':1.64003e-4, 'Ca':2.9756e-4, 'Sc':0.0, 'Ti':2.8645e-4, 'V':9.5492e-07, 'Cr':1.3287e-4, 'Mn':0.0, 'Fe':8.2444e-05, 'Co':0.0, 'Ni':4.3071e-4, 'Cu':2.10858e-4, 'Zn':5.9594e-4, 'Ga':1.9713e-4, 'Ge':5.87597e-4, 'As':0.0, 'Se':4.6268e-4, 'Br':1.56275e-4, 'Kr':2.4849e-4,
             'Rb':1.0969e-4, 'Sr':6.0994e-05, 'Y':0.0, 'Zr':3.42626e-4, 'Nb':0.0, 'Mo':5.9793e-4, 'Tc':0.0, 'Ru':4.0663e-4, 'Rh':0.0, 'Pd':3.0945e-4, 'Ag':8.5796e-05, 'Cd':2.7161e-4, 'In':1.2456e-05, 'Sn':3.34085e-4, 'Sb':6.6075e-05, 'Te':2.8395e-4, 'I':0, 'Xe':2.6779e-4,
             'Cs':0.0, 'Ba':6.2368e-05, 'La':4.7603e-08, 'Ce':2.2495e-05, 'Pr':0.0, 'Nd':2.3159e-4, 'Pm':0.0, 'Sm':3.3472e-4, 'Eu':4.32889e-05, 'Gd':1.27677e-4, 'Tb':0.0, 'Dy':5.20756e-05, 'Ho':0.0, 'Er':7.2459e-05, 'Tm':0.0, 'Yb':8.5449e-05, 'Lu':8.2759e-07, 'Hf':5.2536e-05,
             'Ta':3.80667e-09, 'W':6.9669e-05, 'Re':2.7084e-05,'Os':7.4520e-05, 'Ir':2.5378e-05, 'Pt':3.39199e-05, 'Au':0.0, 'Hg':6.5260e-05, 'Tl':1.99668e-05, 'Pb':1.94476e-05, 'Bi':0.0}

def calc_MFPS(Elem_tabl):
    
    tab_len  = len(Elem_tabl)
    Mass = [pmg.Element[Elem_tabl[i]].atomic_mass for i in np.arange(tab_len)]
    MassSum = np.sum(Mass)
    MFPS = 0.0
    
    for i in np.arange(tab_len):
        MFPS = MFPS + (Mass[i]/MassSum)**2 * MassFluct[Elem_tabl[i]]
    
    MFPS = tab_len * MFPS
    
    return MFPS

if __name__ == '__main__':
   MFPS = calc_MFPS(['Ga','N'])
#   print("%.4e" % MFPS)