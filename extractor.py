# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\smdes\.spyder2\.temp.py
"""
'''
design of countercurrent reactive extractor:
    for removal of mercaptans: ethanethiol, propanethiol and butanethiol
    from a feed of sour kerosene
    using caustic soda solution as extractant
'''

import math
Qk = float(input("Kerosene flow rate in m3/hr: "))
Qw = float(input("Water flow rate in m3/hr: "))

dk = 800.0 # kg/m3 density of kerosene
dw = 1000.0 # density of water

nk = 0.00164 #viscosity
nw = 0.001

ift = 0.052 # interfacial tension of water and kerosene in mN/m

Ce1 = input("concentration of EtSH at inlet ppm: ")
Cp1 = input("concentration of PrSH at inlet ppm: ")
Cb1 = input("concentration of BuSH at inlet ppm: ")
Cm_out = input("required concnetraiton of mercaptans at outlet ppm:")

De = 1.16*(10**-9) # diffusivity of ethanethiol
Dp = 1.01*(10**-9) # diffusivity of propanethiol
Db = 0.911*(10**-9) # diffusivity of butanethiol
'''
The diffusivities were found from Wilke Chang data
'''
dp = 1.15*((ift/((dw-dk)*9.81))**0.5)  
c_bottom = input("weight fraction of base in water: ") # wt.fraction of base NaOH
c= c_bottom-(Ce1+Cb1+Cp1-Cm_out)*(10**-6)
Doh = 5.3*(10**-9) # diffusivity of hydroxyl ions
'''
Distribution Coefficients of mercaptans
'''
Pet = 162.0
Ppr = 50.0
Pbu = 15.0

packing = input("(input one if column is unpacked input 2 if packed ) : ") 
'''
The following if else statements allows users to choose between
getting a packed column or 
an unpacked dispesion column
'''  

if packing == 1:
    ap = 1.0
    trt = 1.0 # totrtuosity
    vdfr = 1.0

else:
    ap = input("packing size: ")
    trt = ap*dp*0.5
    vdfr = input("void factor of packings:")

v =  Qw/(Qk+Qw)  # guess for hold up

if v < 1:
    phi = v
else:
    phi = 0.52

Vk = 0.33*((dw-dk)**0.4)*(dw**-0.6)*(dp**-0.8)*(nk**0.2)*(vdfr/ap)    

Vso = (dw-dk)*dp*dp*9.81/(18*nk)
Vdf = 0.5*Vso
Vcf = 0.2*Vso
cosine = (math.cos(0.25*trt*math.pi))**-2
error = 10

'''
Calculation of Flooding Velocity using Siebert Fair model
'''

while (error > 0.000001):
    Vcf =((0.192*Vso*vdfr)-(cosine)*Vdf)/1.08
    error = (phi - ((Vdf*cosine)/((vdfr*Vso*math.exp(-1.92*phi)) - (Vcf/(1-phi)))))**2
    Vdf = Vdf - (error*0.1)
    phi = phi + error*0.001
    
print "Vcf = ", Vcf
print "Vdf = ", Vdf
print "hold up = ", phi
print "error in hold up = ", error
'''
Calculation of diameter from known flooding velocity
'''
D =(((Qk+Qw)*4)/0.6*(Vcf+Vdf)*math.pi)**0.5 ## taking actual velocity as 0.6 of flooding velocity
Vc = Qk/(math.pi*(D**2)*0.25*3600)
Vd = Qw/(math.pi*(D**2)*0.25*3600)
Re_disp = dw*Vd*dp/nw
Cd = (24/Re_disp)
Vslip = Vk*(1-phi)
Ac = 0.25*(D**2)*math.pi
print "Diameter of Column = ", D
'''
Calculation of mass transfer coefficients
mass transfer cefficient is calculated using the formula:
    Kod = A*A1*A2*A3  where
    A is a constant depending on the mercaptan being trasnferred
    A1 = [(0.95 + c)^4.119] c is mol fraction of base in water phase
    A2 = [(particle size/void fraction)^0.091]
    A3 = [(100*Vd)^0.837] Vd is superficial velocity of disperesed phase
'''

Aet = 4.23    
Apr = 7.756
Abu = 17.65
A1 = (0.95+c)**4.119
A2 = (ap/vdfr)**0.091
A3 = (100*Vd)**0.837
Kod_e = Aet*A1*A2*A3
Kod_b = Abu*A1*A2*A3
Kod_p = Apr*A1*A2*A3

ye = 0
yp = 0
yb = 0

Ce = Ce1
Cp = Cp1
Cb = Cb1
h=0
dh = 0.5
c_err = 1
A = []
B=[]
C=[]
D=[]

'''
using euler's method for dynamic mass balance equations
 to calculate height of column required for 
'''
while abs(c_err) > 0.0001:
    while Ce+Cp+Cb > Cm_out: 
        A.append(Ce)
        B.append(Cp)
        C.append(Cb)
        D.append(h)
        if Ce <= 0:
            E_e = 0
            Ce = 0
            dne = 0
            
        else:
            E_e = 1 + (Doh*(c*40/((c*40 +(1-c)*18)))/(De*Pet*Ce))
            dne = E_e*Kod_e*(Pet*Ce)*Ac*dh*3600*(10**-6)
        
        if Cp <= 0:
            Cp = 0
            dnp = 0
            E_p = 0
        else:
            E_p = 1 + (Doh*c*(c*40/((c*40 +(1-c)*18))))/(Dp*Ppr*Cp)
            dnp = E_p*Kod_p*(Ppr*Cp)*Ac*dh*3600*(10**-6)
    
        if Cb <=0:
            E_b = 0
            Cb = 0
            dne = 0
        else:
            E_b = 1 + (Doh*c*(c*40/((c*40 +(1-c)*18))))/(Db*Pbu*Cb)
            dnb = E_b*Kod_b*(Pbu*Cb)*Ac*dh*3600*(10**-6)
        
        
        if c<=0:
            dne = 0
            dnp = 0
            dnb = 0
            print "Ce = ",Ce
            print "Cb = ",Cb
            print "Cp = ",Cp
            print "Concentration of base is insufficient"
            print h
        Ce = Ce - dne/Qk
        Cp = Cp - dnp/Qk
        Cb = Cb - (dnb/Qk)
        ye = ye+(dne/Qw)
        yp = yp + (dnp/Qw)
        yb = yb + (dnb/Qw)
        c= c+ (dne+dnp+dnb)/(1000000*Qw/40)
        h = h+dh
    c_err = c - c_bottom
    c = c - 0.01*c_err
    
from pylab import plot,show
plot(D,A,'green')
plot(D,B,'blue')
plot(D,C,'red')


print "Concentration of ethanethiol at outlet = ",Ce
print "Concnetration of propanethiol at outlet = ",Cp
print "Concentration of butanethiol at outlet = ",Cb
print "height of column = ",h
print "concentration of base at top of column = ",c
show()