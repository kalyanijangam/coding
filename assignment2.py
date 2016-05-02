# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\admin\.spyder2\.temp.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
import matplotlib.pyplot as plt

m1=1
m2=2
a=327
P=10
U=100
L=10.0
n=10.0


def temp(T,x):
    m1=1
    m2=2
    P=10
    U=100
    cp1=4000+0.1*T[0]+.01*(T[0]**2)
    cp2=3000+0.2*T[1]+.05*(T[1]**2)
    dTa =-(P*U*(T[0]-T[1]))/(m1*cp1)
    dTb =-(P*U*(T[0]-T[1]))/(m2*cp2)
    return [dTa,dTb]
T0=[400,a]
x=np.linspace(0,L,n)
tem=odeint(temp,T0,x)  
print(x)
print(tem[n-1,1])
dx=0.0005
f=1
while (tem[n-1,1]>300.0):
    T0=[400,a-dx*f]
    P=10
    U=100
    tem=odeint(temp,T0,x)
    f=f+1
print(tem)
k1=tem[n-1,0]
k2=tem[0,1]
#plt.plot(x,tem)  
#plt.show()
   
def int1(T):
    cp1=4000+0.1*T+.01*(T**2)
    return cp1
y1,e1=quad(int1,k1,400.0)
#print(y1)
def int2(T):
    cp2=3000+0.2*T+.05*(T**2)
    return cp2
y2,e2=quad(int2,300.0,k2)
#print(y2)
errp=(m1*y1-m2*y2)*100.0/(m1*y1)
print (errp)




 

   



    
    
    
    
