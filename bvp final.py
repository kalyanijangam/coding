# -*- coding: utf-8 -*-
"""
Created on Sun May 01 21:30:10 2016

@author: Archit Datar
"""

''' This is code to solve the 2nd order BVP d2y/dx2= -3*dy/dx+ 2*y with boundary conditions @ x=0, dy/dx= 1
and @ x=2 y= 3
'''

import numpy as np
import scipy 
import scipy.optimize as sco
import scipy.integrate as sci
import matplotlib.pyplot as plt

def func(y,x):
    dydx= np.zeros_like(y)
    dydx[0]= y[1]
    dydx[1]= -3* y[1]+2*y[0]
    return dydx

  
def error(y_init, y_dash_init, y_final):
    y0 =[y_init,y_dash_init]
    y_cal= sci.odeint(func,y0,x)
    error= y_final- y_cal[-1,0]
    return error
    
x= np.linspace(0,2,50)
y_init=0 #initial guess

y_dash_init= 1; y_final=3 #speicifed it question

res= sco.leastsq(error, y_init, args=(y_dash_init,y_final,),full_output=1)
y_init= res[0]
y_error= scipy.sqrt(res[1].diagonal())

#getting everyithing else
y0 =[y_init,y_dash_init]
y_cal= sci.odeint(func,y0,x)

plt.style.use('dark_background')
fig= plt.figure()
ax= fig.add_subplot(111)
ax.plot(x, y_cal)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.title.set_text('BVP')
ax.legend(['y','y_dash'])
