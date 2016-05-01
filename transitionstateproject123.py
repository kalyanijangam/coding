# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:05:44 2016

@author: Hurshvardhai
"""

"""We aim to model the free energy simulation of an SN1 and SN2 reaction
based on breaking and formation of covalent bonds by studying the morse 
potential energy characteristics of the bonds concerned!!"""
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import fsolve
"""We shall first generate the free energy curves for both C-Cl bond and
C-OH bonds with values from literature using classes."""
class morse:
    we=1.0#vibrational const
    wexe=1.0#anharmonicity constant
    c=299792458#speed of light
    h=6.626*(10**(-34))#Js
    Re=1.0#equilibrium bond length
    Na=6.023*(10**23)
    eV=1.602*(10**(-19))#J
    x=1.0
    n=0#vibrational state
    mc=12/Na#g
    mo=35.5
    startr=0.108*(10**(-9))
    endr=1.5*((10**(-9)))
    def morseit(self):
        we=self.we
        wexe=self.wexe
        c=self.c
        h=self.h
        Re=self.Re
        Na=self.Na
        eV=self.eV
        x=self.x
        mc=self.mc
        mo=self.mo
        n=self.n
        startr=self.startr
        endr=self.endr
        mo1=mo/Na
        meu=((mc*mo1)/(mc+mo1))*0.001#in kg
        Do=x*eV#energy at 0 vibration(minimal energy)
        Evib=h*c*((n+0.5)*we-((n+0.5)**2)*wexe)
        De=Do+Evib
        beta=2*3.142*c*we*((meu/(2*De))**0.5)
        r = np.linspace(startr,endr,100)
        print r
        ME=De*((1-(scipy.exp(-beta*(r-Re))))**2)#Morse energy
        ME=ME-De*((1-(scipy.exp(-beta*(endr-Re))))**2)
        Mem=ME*Na/1000#Morse energy mer mole in kJ
        print Mem
        print Do
        print De
        print Re
        print meu
        print beta
        
        plt.plot(r,Mem)
        plt.title("Plot of energy in kJ/mol vs distance in m")
        plt.show()
        return De
"""The class we have created will give us morse characteristics for any bond
given the concerned parameters :)"""
chlorine=morse()
chlorine.mo=35.5
chlorine.Re=1.275*(10**(-10))#m
chlorine.we=2.899*(10**5)#m-1
chlorine.wexe=5205#m-1
chlorine.x=4.436
chlorine.n=0#0 vibrational state while being in the bond
e=chlorine.morseit()
""" Now,
for OH"""
alcohol=morse()
alcohol.mo=17
alcohol.Re=1.012*(10**(-10))#m
alcohol.we=2.899*(10**5)#m-1
alcohol.wexe=5205#m-1
alcohol.x=6.954
alcohol.startr=0.07*(10**(-9))
alcohol.endr=1.2*((10**(-9)))
e=alcohol.morseit()

""" So initially the C is bonded to Cl at a distance equal
to its bond length. OH approaches and puts C-Cl into higher vibrational
states, starting from higher states of its own. After reaction completion,
C-OH is in n=0 and C-Cl is in n=infi but from the modelled eqn, this infinity comes at 
around n=27 by taking derivative and locating the maximum"""

def rfindercl(rguess):
    return(MEcl-Decl*((1-(scipy.exp(-beta*(rguess-chlorine.Re))))**2))
def rfinderoh(rguess):
    return(MEoh-Decl*((1-(scipy.exp(-beta*(rguess-alcohol.Re))))**2))
rc=chlorine.Re 
ra=chlorine.Re+(chlorine.Re-alcohol.Re)     
t=0
dcl=1.275*(10**(-10))  
ncl=0
noh=27
Docl=(chlorine.x)*(chlorine.eV)#energy at 0 vibration(minimal energy)
Dooh=(alcohol.x)*(alcohol.eV)
t=0
while ncl<27:
    ncl=ncl+1
    Evibcl=(chlorine.h)*(chlorine.c)*((ncl+0.5)*(chlorine.we)-((ncl+0.5)**2)*(chlorine.wexe))
    Decl=Docl+Evibcl
    meu=((chlorine.mc*((chlorine.mo)/(chlorine.Na))/(chlorine.mc+((chlorine.mo)/(chlorine.Na)))))*0.001
    beta=2*3.142*(chlorine.c)*(chlorine.we)*((meu/(2*Decl))**0.5)
    MEcl=Decl*((1-(scipy.exp(-beta*(ra-chlorine.Re))))**2)#Morse energy
    ra=chlorine.Re-((scipy.log(1-((MEcl/Decl)**0.5)))/beta)
    MEcl=MEcl*chlorine.Na/1000       
    
    noh=noh-1
    Eviboh=(alcohol.h)*(alcohol.c)*((noh+0.5)*(alcohol.we)-((noh+0.5)**2)*(alcohol.wexe))
    Deoh=Dooh+Eviboh
    meu=((alcohol.mc*((alcohol.mo)/(alcohol.Na))/(alcohol.mc+((alcohol.mo)/(alcohol.Na)))))*0.001
    beta=2*3.142*(chlorine.c)*(chlorine.we)*((meu/(2*Deoh))**0.5)
    MEoh=Deoh*((1-(scipy.exp(-beta*(rc-alcohol.Re))))**2)#Morse energy
    rc=alcohol.Re-((scipy.log(1-((MEoh/Deoh)**0.5)))/beta)
    MEoh=MEoh*chlorine.Na/1000
    
    t=t+1
    ME=MEcl+MEoh
    print(ME,"kJ/mol")
    plt.plot(t,ME,'ro')
    plt.title("Free energy curve for SN2 reaction")

plt.show()

rc=chlorine.Re 
ra=chlorine.Re+(chlorine.Re-alcohol.Re)     
t=0
dcl=1.275*(10**(-10))  
ncl=0
noh=40
Docl=(chlorine.x)*(chlorine.eV)#energy at 0 vibration(minimal energy)
Dooh=(alcohol.x)*(alcohol.eV)
t=0
""" For an SN1 reaction, one bond breakage completes and only then the next 
one begins and thus, we work with two loops"""
while ncl<40:
    ncl=ncl+1
    Evibcl=(chlorine.h)*(chlorine.c)*((ncl+0.5)*(chlorine.we)-((ncl+0.5)**2)*(chlorine.wexe))
    Decl=Docl+Evibcl
    meu=((chlorine.mc*((chlorine.mo)/(chlorine.Na))/(chlorine.mc+((chlorine.mo)/(chlorine.Na)))))*0.001
    beta=2*3.142*(chlorine.c)*(chlorine.we)*((meu/(2*Decl))**0.5)
    MEcl=Decl*((1-(scipy.exp(-beta*(ra-chlorine.Re))))**2)#Morse energy
    ra=chlorine.Re-((scipy.log(1-((MEcl/Decl)**0.5)))/beta)
    MEcl=MEcl*chlorine.Na/1000       
    
    
    
    t=t+1
    t0=t
    ME=MEcl
    print(ME,"kJ/mol")
    plt.plot(t,ME,'ro')
   
"""Now that the intermediate has developed, a resting phase is needed as it is
stable for a very small amount of time"""
s=0
while s<4:
    ME=(t-t0)*(t-t0-4)+MEcl
    plt.plot(t,ME,'ro')
    s=s+1
    t=t+1
while noh>0:
    noh=noh-1
    Eviboh=(alcohol.h)*(alcohol.c)*((noh+0.5)*(alcohol.we)-((noh+0.5)**2)*(alcohol.wexe))
    Deoh=Dooh+Eviboh
    meu=((alcohol.mc*((alcohol.mo)/(alcohol.Na))/(alcohol.mc+((alcohol.mo)/(alcohol.Na)))))*0.001
    beta=2*3.142*(chlorine.c)*(chlorine.we)*((meu/(2*Deoh))**0.5)
    MEoh=Deoh*((1-(scipy.exp(-beta*(rc-alcohol.Re))))**2)#Morse energy
    rc=alcohol.Re-((scipy.log(1-((MEoh/Deoh)**0.5)))/beta)
    MEoh=MEoh*chlorine.Na/1000
    if(noh==39):
        me=MEcl-MEoh
    t=t+1
    ME=MEoh+me
    print(ME,"kJ/mol")
    plt.plot(t,ME,'ro')
    plt.title("Free energy curve for SN1 reaction")
plt.show()

"""Thus as expected SN1 reaction goes via a local minima which is not
observed for an SN2 reaction"""