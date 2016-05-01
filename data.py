# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 11:57:49 2015

@author: Ankit
"""
import scipy
#data file
kla=0.001  #/s
kga=1   #/s
P=101325  #Pa
rho=55.55 #kmol/m3
n=10
L=12
d=4
S=(scipy.pi/4)*(d**2)
Fv_CO2_in=45
Fv_H2O_in=0
Fv_H2S_in=5
Fv_CH4_in=50
Fl_H2O_in=100
Fl_CO2_in=0
Fl_H2S_in=0
Fl_CH4_in=0
Tin_CH4=300+273.16
Tin_H2S=300+273.16
Tin_CO2=300+273.16
Tin_H2O=30+273.16
T_v=300+273.16
T_l=30+273.16