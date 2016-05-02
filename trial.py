# -*- coding: utf-8 -*-
"""
Created on Sun May 01 16:16:18 2016

@author: USER
"""
import scipy
import numpy as np
import scipy.optimize
import pandas as pd
"""T=[0,9,18,26,35,44,53,62,71,80]
A_data=[49.44,41.17,33.54,27.54,20.96,14.93,10.02,7.71,4.8,1.88]
B_data=[0,1.14,0.88,3.44,4.13,4.48,4.75,5.57,4.75,5.8]
p0=np.array([0,0])
def func(T,p):
    return (A_data[0]**(1-p[0])-p[1]*(1-p[0])*T)**(1/(1-p[0]))

def error(p):
    return A_data- func(T,p)
    
res=scipy.optimize.leastsq(error,p0,full_output=1)
p= res[0]
pcov=res[1].diagonal()
A_cal=func(T,p)
print A_cal
"""

e=pd.ExcelFile('C:\\Users\\USER\\Documents\\6th sem\\sim lab python\\work\\ExamProblemData2.1.xlsx')
sheet=e.parse(0)
x=sheet.icol(0).real

y=sheet.icol(1).real

def curve(x,p):
    [A,B,C]=p
    y=A+B*x+(C*(x**2))
    return y
    
def error(p,x,xexp):
    xcalc=curve(x,p)
    err=xexp-xcalc
    return err
    
pguess=[20,30,30]

xlsq=scipy.optimize.leastsq(error,pguess,args=(y,x))

q=xlsq[0]
print q

ycalc=curve(x,q)

print ycalc

error1=ycalc-y
print error1