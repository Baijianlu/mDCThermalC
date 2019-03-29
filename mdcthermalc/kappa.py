# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:11:25 2019

@author: Tao.Fan 
this script used for calculating the integral of kappa
eienstein frequency need to be considered how to calculate, is it equal to Debye frequency?

"""
import os
import numpy as np
import scipy.constants
from scipy.integrate import quad
from mdcthermalc.Get_GVD import Get_GVD
from mdcthermalc.calc_MFPS import calc_MFPS
import pymatgen as pmg

Planck = scipy.constants.hbar
Boltzm = scipy.constants.Boltzmann
atommass = scipy.constants.atomic_mass

def Get_abundancy():
    return 1.256*1e-4                          ###Mg2Si:4.3387*1e-4   LiMg:0.0010  Al3Li:1.256*1e-4

def t_Umklapp(grun,velo,Debye,mass,T):
    return (grun**2 * Boltzm**2 * T**3)/(mass * velo**2 * Debye * Planck) * np.exp(-Debye/(3*T))

def t_Normal(grun,velo,mass,vol,T):
    return (grun**2 * Boltzm**5 * T**5 * vol)/(mass * velo**5 * Planck**4)

def t_Isotope(velo,vol,abund,T):
    return (vol * Boltzm**4 * abund * T**4)/(4 * np.pi * Planck**4 * velo**3)    

def constC(velo):
    return Boltzm**4/(2 * np.pi**2 * Planck**3 * velo)

def get_fun1(x,RT_N,RT_U,RT_ISO):
    return 1/(RT_N * x + RT_U * x**2 + RT_ISO * x**4) * x**4 * np.exp(x)/(np.exp(x)-1)**2

def get_fun2(x,RT_N,RT_U,RT_ISO):
    return RT_N/(RT_N + RT_U * x + RT_ISO * x**3) * x**4 * np.exp(x)/(np.exp(x)-1)**2

def get_fun3(x,RT_N,RT_U,RT_ISO):
    return RT_N * (RT_U + RT_ISO * x**2) /(RT_N + RT_U * x + RT_ISO * x**3) * x**6 * np.exp(x)/(np.exp(x)-1)**2 

def get_fun4(x,RT_N,RT_U,RT_ISO):
    return 1/(RT_N + RT_U + RT_ISO * x**2) * x**2 * np.exp(x)/(np.exp(x)-1)**2

def get_fun5(x,RT_N,RT_U,RT_ISO):
    return RT_N/(RT_N + RT_U + RT_ISO * x**2) * x**4 * np.exp(x)/(np.exp(x)-1)**2

def get_fun6(x,RT_N,RT_U,RT_ISO):
    return RT_N * (RT_U + RT_ISO * x**2)/(RT_N + RT_U + RT_ISO * x**2) * x**6 * np.exp(x)/(np.exp(x)-1)**2

def HeatCapacity(ADebye, ODebye, T, struct):
    N = 1                          # number of primitive cell
    prims = struct.get_primitive_structure()
    Vol = prims.volume * 1e-30             # primitive cell volume  LiMg:44:791067 Mg2Si:64.4719*1e-30
    p = prims.composition.num_atoms                          # atom number in primitive cell
    fun = lambda x: x**4 * np.exp(x)/(np.exp(x)-1)**2
    Cv_aco = 9 * N/Vol * Boltzm * (T/ADebye)**3 * quad(fun,0,ADebye/T)[0]
    Cv_opt = (3*p-3) * N/Vol * Boltzm * (ODebye/T)**2 * np.exp(ODebye/T)/(np.exp(ODebye/T)-1)**2
    return Cv_aco, Cv_opt

def Kappa(filepath,Temp=300.0):
    struct = pmg.Structure.from_file(filepath + '/POSCAR')
    M_avg = 0.0
    for ele in struct.symbol_set:
        M_avg = M_avg + pmg.Element(ele).atomic_mass * struct.composition.get_atomic_fraction(ele)
    
    M_avg = atommass * M_avg             ##kg  Mg2Si:25.5625   LiMg:15.623   Al3Li:21.9714
    V_avg = struct.volume/struct.composition.num_atoms * 1e-30                 ##m^3 Mg2Si:21.49063*1e-30   LiMg:22.3955   Al3Li:16.1203*1e-30
    
    (gruneisen,velocity,DebyeT,freq) = Get_GVD(filepath)
    velocity = velocity * 1e2
    abund = calc_MFPS(list(struct.symbol_set))
    kappa = np.zeros(4)
    avgkappa = np.zeros(len(Temp))
    ADebye = np.sum(DebyeT[0:3]*freq[0:3])/np.sum(freq[0:3])
    ODebye = DebyeT[3]
    
    fp = open('kappa','w')
    fp.write('Temp[K]     Kappa[W*m/K]     R_A     R_O%s' % os.linesep)
    
    for k in np.arange(len(Temp)):
        T = Temp[k]
        for branch in np.arange(4):         # three acoustic branch and one optic branch
            if branch == 0 or branch ==1:
                t_TU = t_Umklapp(gruneisen[branch],velocity[branch],DebyeT[branch],M_avg,T)
                t_TN = t_Normal(gruneisen[branch],velocity[branch],M_avg,V_avg,T)
                t_TISO = t_Isotope(velocity[branch],V_avg,abund,T)
                C_T = constC(velocity[branch])
                IT_1 = C_T * T**3 * quad(get_fun1, 0.0, DebyeT[branch]/T, args=(t_TN,t_TU,t_TISO))[0]
                BettaT_1 = quad(get_fun2, 0.0, DebyeT[branch]/T, args=(t_TN,t_TU,t_TISO))[0]
                BettaT_2 = quad(get_fun3, 0.0, DebyeT[branch]/T, args=(t_TN,t_TU,t_TISO))[0]
                IT_2 = C_T * T**3 * BettaT_1**2/BettaT_2
                kappa[branch] = IT_1 + IT_2
            elif branch == 2:
                t_LU = t_Umklapp(gruneisen[branch],velocity[branch],DebyeT[branch],M_avg,T)
                t_LN = t_Normal(gruneisen[branch],velocity[branch],M_avg,V_avg,T)
                t_LISO = t_Isotope(velocity[branch],V_avg,abund,T)
                C_L = constC(velocity[branch])
                IL_1 = C_L * T**3 * quad(get_fun4, 0.0, DebyeT[branch]/T, args=(t_LN,t_LU,t_LISO))[0]
                BettaL_1 = quad(get_fun5, 0.0, DebyeT[branch]/T, args=(t_LN,t_LU,t_LISO))[0]
                BettaL_2 = quad(get_fun6, 0.0, DebyeT[branch]/T, args=(t_LN,t_LU,t_LISO))[0]
                IL_2 = C_L * T**3 * BettaL_1**2/BettaL_2
                kappa[branch] = IL_1 + IL_2                        
            else:
                t_OU = t_Umklapp(gruneisen[branch],velocity[branch],DebyeT[branch],M_avg,T)
                t_ON = t_Normal(gruneisen[branch],velocity[branch],M_avg,V_avg,T)
                t_OISO = t_Isotope(velocity[branch],V_avg,abund,T)
                C_O = constC(velocity[branch])
                IO_1 = C_O * T**3 * quad(get_fun4, 0.0, DebyeT[branch]/T, args=(t_ON,t_OU,t_OISO))[0]
                BettaO_1 = quad(get_fun5, 0.0, DebyeT[branch]/T, args=(t_ON,t_OU,t_OISO))[0]
                BettaO_2 = quad(get_fun6, 0.0, DebyeT[branch]/T, args=(t_ON,t_OU,t_OISO))[0]
                IO_2 = C_O * T**3 * BettaO_1**2/BettaO_2
                kappa[branch] = IO_1 + IO_2                         
        (Cv_a, Cv_o) = HeatCapacity(ADebye,ODebye,T,struct)
        Rat = Cv_a/(Cv_a+Cv_o)
        avgkappa[k] = Rat * np.average(kappa[0:3]) + (1-Rat) * kappa[3]
        
        fp.write('%-12.1f%-17.3f%-8.3f%-.3f%s' % (T,avgkappa[k],Rat,1-Rat,os.linesep))
    
    fp.close()
    return avgkappa

if __name__ == '__main__':
   avgkappa = Kappa("D:\\nwpuf\\Documents\\mDCThermal\\MgSi",[300.0,400.0,500.0,600.0,700.0,800.0,900.0])    

