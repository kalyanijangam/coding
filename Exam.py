# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 11:55:44 2015

@author: Ankit
"""

import scipy 
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy import integrate
import data
from Compounds import CO2,H2S,CH4,H2O
class Exam():
    def __init__(self):
        #initialize all ur inputs
        self.Fv_CO2_in=data.Fv_CO2_in
        self.Fv_H2O_in=data.Fv_H2O_in
        self.Fv_H2S_in=data.Fv_H2S_in
        self.Fv_CH4_in=data.Fv_CH4_in
        self.Fl_CO2_in=data.Fl_CO2_in
        self.Fl_H2O_in=data.Fl_H2O_in
        self.Fl_H2S_in=data.Fl_H2S_in
        self.Fl_CH4_in=data.Fl_CH4_in
        self.Tin_CO2=data.Tin_CO2
        self.Tin_H2O=data.Tin_H2O
        self.Tin_H2S=data.Tin_H2S
        self.Tin_CH4=data.Tin_CH4
        self.T_v=data.T_v
        self.T_l=data.T_l
        self.kla=data.kla
        self.kga=data.kga
        self.P=data.P
        self.n=data.n
        self.L=data.L
        self.d=data.d
        self.S=data.S
        self.rho=data.rho
     
        
    def set_grid(self):
        #Guess values array for flow rate
        self.Fl_CO2=scipy.ones(self.n)*self.Fv_CO2_in
        self.Fl_H2S=scipy.ones(self.n)*self.Fv_H2S_in
        self.Fl_CH4=scipy.ones(self.n)*self.Fv_CH4_in
        self.Fl_H2O=scipy.ones(self.n)*self.Fl_H2O_in
        self.Fv_CO2=scipy.ones(self.n)*self.Fv_CO2_in
        self.Fv_H2O=scipy.ones(self.n)*self.Fl_H2O_in
        self.Fv_H2S=scipy.ones(self.n)*self.Fv_H2S_in
        self.Fv_CH4=scipy.ones(self.n)*self.Fv_CH4_in
        
        #guess value for enthalpy
        self.Hin_CO2=self.Fv_CO2_in*(scipy.integrate.quad(CO2.CpG,298.15,self.Tin_CO2)[0]+CO2.Hf)
        self.Hin_CH4=self.Fv_CH4_in*(scipy.integrate.quad(CH4.CpG,298.15,self.Tin_CH4)[0]+CH4.Hf)
        self.Hin_H2S=self.Fv_H2S_in*(scipy.integrate.quad(H2S.CpG,298.15,self.Tin_H2S)[0]+H2S.Hf)
        self.Hin_H2O=self.Fv_H2O_in*(scipy.integrate.quad(H2O.CpG,298.15,self.Tin_H2O)[0]+H2O.Hf-H2O.Hv(H2O.Tb))
        
        """self.H_CH4=scipy.ones(self.n)*self.Hin_CH4
        self.H_CO2=scipy.ones(self.n)*self.Hin_CO2
        self.H_H2S=scipy.ones(self.n)*self.Hin_H2S
        self.H_H2O=scipy.ones(self.n)*self.Hin_H2O"""
        
        self.H_v=scipy.ones(self.n)*self.Hin_CH4
        self.H_l=scipy.ones(self.n)*self.Hin_CH4
        
        #guess values for Temp
        self.T=scipy.ones(self.n)*self.T_v
        self.dz=self.L/float(self.n-1)
        
        
        #concatenate alll ur guesses
        self.F_guess=scipy.concatenate((self.Fl_H2O,self.Fl_H2S,self.Fl_CO2,self.Fl_CH4,self.Fv_CH4,self.Fv_H2O,self.Fv_H2S,self.Fv_CO2,self.T))
        """self.H_guess=scipy.concatenate((self.H_H2O,self.H_H2S,self.H_CO2,self.H_CH4))
        self.guess=scipy.concatenate((self.F_guess,self.H_guess))"""
        
    def solve(self):
        F_guess=self.F_guess
        self.F_H_soln=scipy.optimize.leastsq(residuals,F_guess,args=(self))[0]
        
        self.Fl_H2O=self.F_H_soln[:self.n]
        self.Fl_H2S=self.F_H_soln[self.n:2*self.n]
        self.Fl_CO2=self.F_H_soln[2*self.n:3*self.n]
        self.Fl_CH4=self.F_H_soln[3*self.n:4*self.n]
        
        self.Fv_CH4=self.F_H_soln[4*self.n:5*self.n]
        self.Fv_H2O=self.F_H_soln[5*self.n:6*self.n]
        self.Fv_H2S=self.F_H_soln[6*self.n:7*self.n]
        self.Fv_CO2=self.F_H_soln[7*self.n:8*self.n]
        
        """self.H_H2O=self.F_H_soln[8*self.n:9*self.n]
        self.H_H2S=self.F_H_soln[9*self.n:10*self.n]
        self.H_CO2=self.F_H_soln[10*self.n:11*self.n]
        self.H_CH4=self.F_H_soln[11*self.n:12*self.n]"""
        
        self.T=self.F_H_soln[8*self.n:9*self.n]
        print "Fv_H2S"
        print self.Fv_H2S
        print "Fl_H2S"
        print self.Fl_H2S
        print "Fv_CO2"
        print self.Fv_CO2
        print "Fl_CO2"
        print self.Fl_CO2
        print "Fv_CH4"
        print self.Fv_CH4
        print "Fl_CH4"
        print self.Fl_CH4
        print "Fl_H2O"
        print self.Fl_H2O
        print "Fv_H2O"
        print self.Fv_H2O
        print "T"
        print self.T
        
def residuals(F,obj):
    kla=obj.kla
    kga=obj.kga
    P=obj.P
    n=obj.n
    L=obj.L
    dz=obj.dz
    d=obj.d
    S=obj.S
    rho=obj.rho
    T_v=obj.T_v
    Fv_CO2_in=obj.Fv_CO2_in
    Fv_H2S_in=obj.Fv_H2S_in
    Fv_CH4_in=obj.Fv_CH4_in
    Fv_H2O_in=obj.Fv_H2O_in
    Fl_CO2_in=obj.Fl_CO2_in
    Fl_H2S_in=obj.Fl_H2S_in
    Fl_CH4_in=obj.Fl_CH4_in
    Fl_H2O_in=obj.Fl_H2O_in
    Tin_CO2=obj.Tin_CO2
    Tin_H2O=obj.Tin_H2O
    Tin_H2S=obj.Tin_H2S
    Tin_CH4=obj.Tin_CH4
    #doing conc cal.
    Fl_H2O=F[:n]
    Fl_H2S=F[n:2*n]
    Fl_CO2=F[2*n:3*n]
    Fl_CH4=F[3*n:4*n]
    
    Fv_CH4=F[4*n:5*n]
    Fv_H2O=F[5*n:6*n]
    Fv_H2S=F[6*n:7*n]
    Fv_CO2=F[7*n:8*n]
    
    #enthalpy terms
    """H_H2O=F[8*n:9*n]
    H_H2S=F[9*n:10*n]
    H_CO2=F[10*n:11*n]
    H_CH4=F[11*n:12*n]"""
    #Temp
    T=F[8*n:9*n]
    
    Fv=scipy.ones(n) #kmol/s
    Fl=scipy.ones(n)*100
    Fv_H2S[0]=5*Fv[0]
    Fv_CH4[0]=50*Fv[0]
    Fv_CO2[0]=45*Fv[0]
    Fv_H2O[0]=0
    Fl_H2S[-1]=0
    Fl_CO2[-1]=0
    Fl_CH4[-1]=0
    Fl_H2O[-1]=100
    
    #for storing errors
    err_H2S_v=scipy.ones(n)
    err_H2S_l=scipy.ones(n)
    err_CO2_v=scipy.ones(n)
    err_CO2_l=scipy.ones(n)
    err_CH4_v=scipy.ones(n)
    err_CH4_l=scipy.ones(n)
    err_H2O_v=scipy.ones(n)
    err_H2O_l=scipy.ones(n)
    
    #Henrys constants
    T=scipy.ones(n)*T_v
    T[0]=300+273.16   
    T[-1]=30+273.16
    
    #forward difference    
    Fv[0]=Fv_H2S[0]+Fv_CO2[0]+Fv_CH4[0]+Fv_H2O[0]
    Fl[0]=Fl_H2S[0]+Fl_CO2[0]+Fl_CH4[0]+Fl_H2O[0]
    err_H2S_v[0]=((1/kga+H2S.Hen(T[0])/kla)/S)*((Fv_H2S[1]-Fv_H2S[0])/dz)+(P*Fv_H2S[0]/Fv[0])-(H2S.Hen(T[0])*Fl_H2S[0]*rho/Fl[0])
    err_CO2_v[0]=((1/kga+CO2.Hen(T[0])/kla)/S)*((Fv_CO2[1]-Fv_CO2[0])/dz)+(P*Fv_CO2[0]/Fv[0])-(CO2.Hen(T[0])*Fl_CO2[0]*rho/Fl[0])
    err_CH4_v[0]=((1/kga+CH4.Hen(T[0])/kla)/S)*((Fv_CH4[1]-Fv_CH4[0])/dz)+(P*Fv_CH4[0]/Fv[0])-(CH4.Hen(T[0])*Fl_CH4[0]*rho/Fl[0])
    err_H2S_l[0]=(Fv_H2S[1]-Fv_H2S[0])/dz-(Fl_H2S[1]-Fl_H2S[0])/dz
    err_CO2_l[0]=(Fv_CO2[1]-Fv_CO2[0])/dz-(Fl_CO2[1]-Fl_CO2[0])/dz
    err_CH4_l[0]=(Fv_CH4[1]-Fv_CH4[0])/dz-(Fl_CH4[1]-Fl_CH4[0])/dz
    err_H2O_v[0]=((1/kga+H2O.Hen/kla)/S)*((Fv_H2O[1]-Fv_H2O[0])/dz)+((P*Fv_H2O[0]/Fv[0])-H2O.Psat(T[0]))
    err_H2O_l[0]=((1/kga+H2O.Hen/kla)/S)*((Fl_H2O[1]-Fl_H2O[0])/dz)+((P*Fl_H2O[0]/Fl[0])-H2O.Psat(T[0]))
    
    #backward difference
    Fv[-1]=Fv_H2S[-1]+Fv_CO2[-1]+Fv_CH4[-1]+Fv_H2O[-1]
    Fl[-1]=Fl_H2S[-1]+Fl_CO2[-1]+Fl_CH4[-1]+Fl_H2O[-1]
    err_H2S_v[-1]=((1/kga+H2S.Hen(T[-1])/kla)/S)*((Fv_H2S[-1]-Fv_H2S[-2])/dz)+(P*Fv_H2S[-1]/Fv[-1])-(H2S.Hen(T[-1])*Fl_H2S[-1]*rho/Fl[-1])
    err_CO2_v[-1]=((1/kga+CO2.Hen(T[-1])/kla)/S)*((Fv_CO2[-1]-Fv_CO2[-2])/dz)+(P*Fv_CO2[-1]/Fv[-1])-(CO2.Hen(T[-1])*Fl_CO2[-1]*rho/Fl[-1])
    err_CH4_v[-1]=((1/kga+CH4.Hen(T[-1])/kla)/S)*((Fv_CH4[-1]-Fv_CH4[-2])/dz)+(P*Fv_CH4[-1]/Fv[-1])-(CH4.Hen(T[-1])*Fl_CH4[-1]*rho/Fl[-1])
    err_H2S_l[-1]=(Fv_H2S[-1]-Fv_H2S[-2])/dz-(Fl_H2S[-1]-Fl_H2S[-2])/dz
    err_CO2_l[-1]=(Fv_CO2[-1]-Fv_CO2[-2])/dz-(Fl_CO2[-1]-Fl_CO2[-2])/dz
    err_CH4_l[-1]=(Fv_CH4[-1]-Fv_CH4[-2])/dz-(Fl_CH4[-1]-Fl_CH4[-2])/dz
    err_H2O_v[-1]=((1/kga+H2O.Hen/kla)/S)*((Fv_H2O[-2]-Fv_H2O[-1])/dz)+((P*Fv_H2O[-1]/Fv[-1])-H2O.Psat(T[-1]))
    err_H2O_l[-1]=((1/kga+H2O.Hen/kla)/S)*((Fl_H2O[-2]-Fl_H2O[-1])/dz)+((P*Fl_H2O[-1]/Fl[-1])-H2O.Psat(T[-1]))
    
    #center difference
    for i in range(1,n-1):
        Fv[i]=Fv_H2S[i]+Fv_CO2[i]+Fv_CH4[i]+Fv_H2O[i]
        Fl[i]=Fl_H2S[i]+Fl_CO2[i]+Fl_CH4[i]+Fl_H2O[i]
        err_H2S_v[i]=((1/kga+H2S.Hen(T[i])/kla)/S)*((Fv_H2S[i+1]-Fv_H2S[i-1])/(2*dz))+(P*Fv_H2S[i]/Fv[i])-(H2S.Hen(T[i])*Fl_H2S[i]*rho/Fl[i])
        err_CO2_v[i]=((1/kga+CO2.Hen(T[i])/kla)/S)*((Fv_CO2[i+1]-Fv_CO2[i-1])/(2*dz))+(P*Fv_CO2[i]/Fv[i])-(CO2.Hen(T[i])*Fl_CO2[i]*rho/Fl[i])
        err_CH4_v[i]=((1/kga+CH4.Hen(T[i])/kla)/S)*((Fv_CH4[i+1]-Fv_CH4[i-1])/(2*dz))+(P*Fv_CH4[i]/Fv[i])-(CH4.Hen(T[i])*Fl_CH4[i]*rho/Fl[i])
        err_H2S_l[i]=(Fv_H2S[i+1]-Fv_H2S[i-1])/(2*dz)-(Fl_H2S[i+1]-Fl_H2S[i-1])/(2*dz)
        err_CO2_l[i]=(Fv_CO2[i+1]-Fv_CO2[i-1])/(2*dz)-(Fl_CO2[i+1]-Fl_CO2[i-1])/(2*dz)
        err_CH4_l[i]=(Fv_CH4[i+1]-Fv_CH4[i-1])/(2*dz)-(Fl_CH4[i+1]-Fl_CH4[i-1])/(2*dz)
        err_H2O_v[i]=((1/kga+H2O.Hen/kla)/S)*((Fv_H2O[i+1]-Fv_H2O[i-1])/(2*dz))+(P*Fv_H2O[i]/Fv[i])-(H2O.Psat(T[i]))
        err_H2O_l[i]=((1/kga+H2O.Hen/kla)/S)*((Fl_H2O[i+1]-Fl_H2O[i-1])/(2*dz))+(P*Fl_H2O[i]/Fv[i])-(H2O.Psat(T[i]))
    
    """Molar conc ends, enthslpy starts to get Temp. whooooohooooooooooo"""
    H_v=scipy.ones(n)
    H_l=scipy.ones(n)
    H_v[0]=Fv_CO2[0]*((scipy.integrate.quad(CO2.CpG,298.15,Tin_CO2)[0])+CO2.Hf)+Fv_CH4[0]*(scipy.integrate.quad(CH4.CpG,298.15,Tin_CH4)[0]+CH4.Hf)+Fv_H2S[0]*(scipy.integrate.quad(H2S.CpG,298.15,Tin_H2S)[0]+H2S.Hf)+Fv_H2O[0]*(scipy.integrate.quad(H2O.CpG,298.15,T[0])[0]+H2O.Hf)
    H_l[0]=Fl_CO2[0]*((scipy.integrate.quad(CO2.CpG,298.15,T[0])[0])+CO2.Hf)-(CO2.Habs)+Fl_CH4[0]*(scipy.integrate.quad(CH4.CpG,298.15,T[0])[0]+CH4.Hf)-(CH4.Habs)+Fl_H2S[0]*(scipy.integrate.quad(H2S.CpG,298.15,T[0])[0]+H2S.Hf)-(H2S.Habs)+Fl_H2O[0]*(scipy.integrate.quad(H2O.CpG,298.15,T[0])[0]+H2O.Hf)-(H2O.Hv(H2O.Tb))
    H_v[-1]=Fv_CO2[-1]*((scipy.integrate.quad(CO2.CpG,298.15,T[-1])[0])+CO2.Hf)+Fv_CH4[-1]*(scipy.integrate.quad(CH4.CpG,298.15,T[-1])[0]+CH4.Hf)+Fv_H2S[-1]*(scipy.integrate.quad(H2S.CpG,298.15,T[-1])[0]+H2S.Hf)+Fv_H2O[-1]*(scipy.integrate.quad(H2O.CpG,298.15,T[0])[0]+H2O.Hf)
    H_l[-1]=Fl_CO2[-1]*((scipy.integrate.quad(CO2.CpG,298.15,T[-1])[0])+CO2.Hf)-(CO2.Habs)+Fl_CH4[-1]*(scipy.integrate.quad(CH4.CpG,298.15,T[-1])[0]+CH4.Hf)-(CH4.Habs)+Fl_H2S[-1]*(scipy.integrate.quad(H2S.CpG,298.15,T[-1])[0]+H2S.Hf)-(H2S.Habs)+Fl_H2O[-1]*(scipy.integrate.quad(H2O.CpG,298.15,Tin_H2O)[0]+H2O.Hf)-(H2O.Hv(H2O.Tb))
    for i in range(1,n-1):
        H_v[i]=Fv_CO2[i]*(scipy.integrate.quad(CO2.CpG,298.15,T[i])[0]+CO2.Hf)+Fv_CH4[i]*(scipy.integrate.quad(CH4.CpG,298.15,T[i])[0]+CH4.Hf)+Fv_H2S[i]*(scipy.integrate.quad(H2S.CpG,298.15,T[i])[0]+H2S.Hf)+Fv_H2O[i]*(scipy.integrate.quad(H2O.CpG,298.15,T[i])[0]+H2O.Hf)
        H_l[i]=Fl_CO2[i]*(scipy.integrate.quad(CO2.CpG,298.15,T[i])[0]+CO2.Hf)-(CO2.Habs)+Fl_CH4[i]*(scipy.integrate.quad(CH4.CpG,298.15,T[i])[0]+CH4.Hf)-(CH4.Habs)+Fl_H2S[i]*(scipy.integrate.quad(H2S.CpG,298.15,T[i])[0]+H2S.Hf)-(H2S.Habs)+Fl_H2O[i]*(scipy.integrate.quad(H2O.CpG,298.15,T[i])[0]+H2O.Hf)-(H2O.Hv(H2O.Tb))
    
    err_H=scipy.ones(n) 
    err_H=(scipy.gradient(H_v)/dz)-(scipy.gradient(H_l)/dz)
    
    err=scipy.concatenate((err_H2S_v,err_H2S_l,err_CO2_v,err_CO2_l,err_CH4_v,err_CH4_l,err_H2O_v,err_H2O_l,err_H))
    #print err
    return err
    x=[1,2,3,4,5,6,7,8,9,10]
    plt.plot(x,Fv_CO2)
    plt.show()   