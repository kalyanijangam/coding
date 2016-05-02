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
class temperature:
    m1=1
    m2=2
    
    P=10
    U=100
    L=10.0
    n=10.0
    T1=400
    T20=350
    T2=335
    
    def errp(self):
        P=self.P
        U=self.U
        T20=self.T20
        L=self.L
        n=self.n
        b=L/(n-1)
        m1=self.m1
        m2=self.m2
        T1=self.T1
        T2=self.T2
        f=0
        dx=.001
        while (T20>300.0):
            T20=T2-dx*f
            k=T20
            T10=T1
            x=1
            #print T20
            #print T10
            #print b
            while (x<n):
                #cp1=(4000+.1*T10+.01*T10**2)
                #cp2=(3000+.2*T20+.05*T20**2)
                h1=(P*U*(T20-T10)/(m1*(4000+.1*T10+.01*T10**2)))
                h2=(P*U*(T20-(T10+(.5*b*h1)))/(m1*(4000+.1*(T10+(.5*b*h1))+.01*(T10+(.5*b*h1))**2)))
                h3=(P*U*(T20-(T10+(.5*b*h2)))/(m1*(4000+.1*(T10+(.5*b*h2))+.01*(T10+(.5*b*h2))**2)))
                h4=(P*U*(T20-(T10+(b*h3)))/(m1*(4000+.1*(T10+(b*h3))+.01*(T10+(b*h3))**2)))
                
                g1=(P*U*(T20-T10)/(m2*(3000+.2*T20+.05*T20**2)))
                g2=(P*U*((T20+(.5*b*g1))-T10)/(m2*(3000+.2*(T20+(.5*b*g1))+.05*(T20+(.5*b*g1))**2)))
                g3=(P*U*((T20+(.5*b*g2))-T10)/(m2*(3000+.2*(T20+(.5*b*g2))+.05*(T20+(.5*b*g2))**2)))
                g4=(P*U*((T20+(b*g3))-T10)/(m2*(3000+.2*(T20+(b*g3))+.05*(T20+(b*g3))**2)))
                
                
                
                Ta=T10+(b/6)*(h1+2*h2+2*h3+h4)
                Tb=T20+(b/6)*(g1+2*g2+2*g3+g4)
                T10=Ta
                T20=Tb
                #print T10,T20
                x+=1
                #print x
            f=f+1
        #print T10
        print T20
        print k 
        #print(tem)
        k1=T10
        k2=k
#plt.plot(x,tem)  
#plt.show()
   
        def int1(T):
            cp1=4000+0.1*T+.01*T**2
            return cp1
        y1,e1=quad(int1,k1,400.0)
        #print(y1)
        def int2(T):
            cp2=3000+0.2*T+.05*(T**2)
            return cp2
        y2,e2=quad(int2,300.0,k2)
        #print(y2)
        err=(m1*y1-m2*y2)*100.0/(m1*y1)
        return err
    








 

   



    
    
    
    
