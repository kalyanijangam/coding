# -*- coding: utf-8 -*-
"""
Spyder Editor


"""
import scipy
import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
'''continuous spray dryer
assumptions: 
1.diameter of particle remains same, water loss from crust similar to drying of polymer beads
2.temperature of solid and water in slurry at any point x remains same Ts=Tw in slurry phase 
3. only water transfer occurs in between two phases and solid is not getting carried with air 
or vice versa there is no air entrainment in slurry '''

# details of problem
l=0.1               #kg/s
Wi=40               #intial water content in percentage
w=0.01*Wi*l         #kg/s
Wn=5                #final water content in percentage
s=l*(1-(0.01*Wi))
Xi=Wi/(100-Wi)      #moisture content on dry basis
Xn=Wn/(100-Wn)
wmol0=(Xi*s)/(18e-03) #initial moles entering in dryer
wmoln=(Xn*s)/(18e-03) #final moles leaving the dryer
Ta=150+273.16         #K  input air temp
air=0.5               #kg/sec
amol=0.5e03/28.97     #gmol/s moles of air
Yi=0
Yn=(l*(Xi-Xn))/air    #overall mass balance
Tn1=331.1             #temp at inlet of slurry
'''data of solid, water and air'''
Cpw=4.186               #kJ/kg K 
Cpa=1.013               #kJ/kg K
Cps=0.96                #kJ/kg K
#Ky=150                  #kg/sq.m hr
dp=300                  #micromtr diameter of polymer bead
Rhop=995                #kg/m3 density of polymer
RhoA=0.946              #kg/m3 density of air in given temp range
'''assumption'''
H=10                    #m height of spray dryer
D0=2                    #m initial assumption for finding the diameter of spray dryer
lmb=2258                #kJ/kg-1 latenet heat of vapourization
MuA=1.983e-05
p=1.01325               #bar
'''log10(P) = A − (B / (T + C))
    P = vapor pressure (bar)
    T = temperature (K)'''
fla=air/RhoA             #air flowrate
Va=fla/(3.14*D0*D0)      #m/s superficial velocity of air
ut=1.75*(((9.81*dp*1e-06*(Rhop-RhoA))/RhoA)**0.5)  #terminal settling velocity of dried particles
Dv=(1.17564e-09*(Tn1**1.75)*1.01325)/p

#gas side mass transfer coeff of water in air 
D=D0
P=3.14*D
A=3.14*D*D
Rep=(RhoA*dp*1e-06*Va)/MuA
Sc=MuA/(RhoA*Dv)
Sh=2+0.6*(Rep**0.5)*(Sc**0.33)
ky=(18*RhoA*Dv*Sh)/(28.87*p*dp*1e-06) #gas side mass transfer coeff of water in air 
U=100
#defining array k=[Tair,Tslurry,CwAIR,CwSLURRY] CwAIR CwSLURRY expressed in moles
def dry(k,x):
     
    Psat=10**((3.55959)-(643.748/(k[0]-198.043)))
    ky=150
    a=1
    dTa =-(P*U*(k[0]-k[1]))/(air*Cpa)
    dTw =-(P*U*(k[0]-k[1]))/(k[3]*Cpw)-(ky*a*A*(Psat-(k[2]/(k[2]+amol)*p))*lmb)
    dMa = ky*a*A*(Psat-(k[2]/(k[2]+amol)*p))
    dMw = -dMa
    #    dTs=-(P*U*(T[0]-T[1]))/(s*Cps)
    return [dTa,dTw,dMa,dMw]
e1=1
e2=1

while (Va>=ut):
    
    while(abs(e1)>0.1,abs(e2)>0.1):
        T00=Ta
        Tn1=331.15 #K
        param=[T00,Tn1,0,wmol0]
        n=10
        z=np.linspace(0,H,n)
        y=odeint(dry,param,z,full_output=1)
        e1=(y[9,0]-333)
        e2=(y[0,3]-wmoln)
        D0=D0+e1*0.1
        return D0
    Va=fla/(3.14*D0*D0)
    D=D0
    P=3.14*D
    A=3.14*D*D
    Rep=(RhoA*dp*1e-06*Va)/MuA
    Sc=MuA/(RhoA*Dv)
    Sh=2+0.6*(Rep**0.5)*(Sc**0.33)
    ky=(18*RhoA*Dv*Sh)/(28.87*p*dp*1e-06)
    return Va,D0
Ds=D0+D0*0.2 #oversizing 20%
print 'diameter of spray dryer=', Ds





