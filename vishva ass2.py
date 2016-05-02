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
class absorption:
    CO2=1
    CH4=2
    H2O=3
    
    y1=.5
    y2=.5
    y3=0
    G0=100 #mol/hr
    G10=y1*G0 #mol/hr
    G20=y2*G0#mol/hr
    G30=y3*G0 #mol/hr
    
    x3=1
    x1=0
    x2=0
    LH=150 #mol/hr
    L1H=x1*LH #mol/hr
    L2H=x2*LH #mol/hr
    L3H=x3*LH #mol/hr
    
    
    
    
    H1= 16000 #atm (https://www.google.co.in/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0ahUKEwi0stu66JfKAhWKiCwKHQRHDAUQFggnMAE&url=http%3A%2F%2Fwebbook.nist.gov%2Fcgi%2Fcbook.cgi%3FID%3DC124389%26Mask%3D10&usg=AFQjCNHl7I5pBiB1hL8K5QeiHfQXovNARA)
    H2= 37037 #bar (https://www.google.co.in/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwiNn9qb6ZfKAhWDQyYKHWb2CQIQFggbMAA&url=http%3A%2F%2Fwebbook.nist.gov%2Fcgi%2Fcbook.cgi%3FID%3DC74828%26Mask%3D10&usg=AFQjCNGoiCVKmUIOo3ApMMpjnt3W8dpVLQ)
    Psat= .003169 #MPa at 25degree celcius (Perry's handbook)
    #Kl1=.163/3600 #m/hr
    #Kl2=.2/3600 #m/hr
    #kg3=1.22/3600 #m/hr
    Kl1=.01
    Kl2=.005
    kg3=2.452*(10**-6)
    H=5 #m
    A=.2 #m2
    a=1 #m-1
    
    
    
    P=1 #atm
    T=25 #degreecelcius
    
    def mol(self):
        G0=self.G0   
        G10=self.G10
        G20=self.G20
        G30=self.G30
        y1=self.y1
        y2=self.y2
        y3=self.y3
        
        LH=self.LH
        L1H=self.L1H
        L2H=self.L2H
        L3H=self.L3H
        
        x1=self.x1
        x2=self.x2
        x3=self.x3
        H1=self.H1
        H2=self.H2
        Psat=self.Psat
        Kl1=self.Kl1
        Kl2=self.Kl2
        kg3=self.kg3
        H=self.H
        A=self.A
        a=self.a
        L1=G10
        L2=G20
        L3=L3H-.1*L3H
        P=self.P
        
        #unit conversion
        H1=H1*18*101325/1000 #Pa*m3/mole        
        H2=H2*18*100000/1000 #Pa*m3/mole
        Psat=Psat*1000000 #Pa
        P=P*101325 #Pa
        #Kl1
        #Kl2
        #kg3
            
        def diff(G,z):
            VL=G[5]*.018/1000
            dG1=Kl1*a*A*((((G[0]*P)/(G[0]+G[1]+G[2]))/H1)-(G[3]/VL))
            dG2=Kl2*a*A*((((G[1]*P)/(G[0]+G[1]+G[2]))/H2)-(G[4]/VL))
            dG3=kg3*a*A*(Psat-((G[2]/(G[2]+G[1]+G[0]))*P))
            dL1=dG1
            dL2=dG2
            dL3=dG3
            return [dG1,dG2,dG3,dL1,dL2,dL3]
        Gi=[G10,G20,G30,L1,L2,L3]
        z=np.linspace(0,H,10)
        y=odeint(diff,Gi,z)
        #print z
        #print y
        dx1=0.1
        err1=1
        err=1
        err2=1
    
        while (abs(err)>.01) or (abs(err1)>.01) or (abs(err2>.01)):
            Gi=[G10,G20,G30,L1,L2,L3] 
            a1=L1
            a2=L2
            a3=L3
            y=odeint(diff,Gi,z)
            err=y[9,5]-L3H
            err1=y[9,4]-L2H
            err2=y[9,3]-L1H
            print err,err1,err2    
            r=max(abs(err),abs(err1),abs(err2))
            
            if r==abs(err):
                L3=a3+dx1*r
            elif r==abs(err1):
                L2=a2-dx1*r
            elif r==abs(err2):
                L1=a1-dx1*r
                
        y1=y[:,0]/(y[:,0]+y[:,1]+y[:,2])
        y2=y[:,1]/(y[:,0]+y[:,1]+y[:,2])
        y1H=(y[9,0])/(y[9,0]+y[9,1]+y[9,2])
        y2H=(y[9,1])/(y[9,0]+y[9,1]+y[9,2])
        y3H=(y[9,2])/(y[9,0]+y[9,1]+y[9,2])
        x10=(y[0,3])/(y[0,3]+y[0,4]+y[0,5])
        x20=(y[0,4])/(y[0,3]+y[0,4]+y[0,5])
        x30=(y[0,5])/(y[0,3]+y[0,4]+y[0,5])
        t1=y[0,3]+y[9,0]
        t2=y[0,4]+y[9,1]
        t3=y[0,5]+y[9,2]
        plt.plot(z,y2,'go')
        plt.plot(z,y1,'ro')
        plt.show()
        print y,"\n",'y1H=',y1H,"\n",'y2H=',y2H,"\n",'y3H=',y3H,"\n",'x10=',x10,"\n",'x20=',x20,"\n",'x30=',x30
           
            
            
            
        

