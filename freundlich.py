# -*- coding: utf-8 -*-
"""
Created on Wed Feb 03 11:49:25 2016

@author: Hurshvardhai
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 00:02:55 2016

@author: Hurshvardhai
"""


import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import statsmodels.stats.stattools as stools
import win32com.client
xl= win32com.client.gencache.EnsureDispatch('Excel.Application')
wb=xl.Workbooks('adsorb.xlsx')
sheet=wb.Sheets('Langmuir')
   
def getdata(sheet, Range):
    data= sheet.Range(Range).Value
    data=scipy.array(data)
    data=data.reshape((1,len(data)))[0]
    return data
ydata=getdata(sheet,"E20:E28")
errdata=getdata(sheet,"C20:C28")
def assign(x, p):
    A,B,C = p
    return (A*(x**2)+B*X+C)
def residuals(p, y, x):
    A,B,C = p
    err = (((y-(A*(x**2)+B*x+C))**2)**0.5)
    return err

x = np.linspace(-1.10,-1.40,9)
#x=[11.3,11.9,12.5,13.1,13.7,14.3,14.9,15.5,16.1]
A=0.7216
B=-1.8472
yreal = (A*x+B)
guessAB = [10, 590]
fitting = leastsq(residuals, guessAB, args=(ydata, x))
popt=fitting[0]
pcov=fitting[1]
print("*******************************")
print("Optimized A and B in AX+B")
print popt[0],popt[1],popt[2]
print("A and B in AX+B FROM RESULTS")
print A,B,C


print("*******************************")
print("FREUNDLICH ISOTHERM")
print("(x/m)=k*C^n")
print('Freundlich K=',A)
print('Freunlich n=',B)
print("*******************************")
yfit=popt[0]*x*x+popt[1]*x+popt[2]

def error_fit(Xdata,popt,pcov):
    Y=popt[0]*x*x+popt[1]*x+popt[2]
    dY=[]
    for i in xrange(len(popt)):
        p=popt[i]
        dp=abs(p)/1e6+1e-20
        popt[i]+=dp
        Yi=popt[0]*x+popt[1]
        dy=(Yi-Y)/dp
        dY.append(dy)
        popt[i]-=dp
        dY=scipy.array(dY)
        A=scipy.dot(dY.T,pcov)
        B=scipy.dot(A,dY)
        sigma2=B.diagonal()
        mean_sigma2=scipy.mean(sigma2)
        M=len(Xdata)
        N=len(popt)
        avg_stddev_data=scipy.sqrt(M*mean_sigma2/N)
        sigma=scipy.sqrt(sigma2)
        return sigma
sig=error_fit(x,popt,pcov)
print("******************")
print("SIGMA VALUES")
print(0.1*sig)
print("******************")
plt.plot(x,assign(x,fitting[0]),'r')
plt.plot(x,ydata,'bo')
y11=popt[0]*x+popt[1]+0.1*sig
plt.plot(x,y11,'g')
y22=popt[0]*x+popt[1]-0.1*sig
plt.plot(x,y22,'g')
plt.errorbar(x,ydata,yerr=errdata,ls="none")
M=len(ydata)
N=len(popt)

yavg=scipy.mean(ydata)


squares=(yfit-yavg)
squaresT=(ydata-yavg)
residuals=(yfit-ydata)

SSM=sum(squares**2)
SSE=sum(residuals**2)
SST=sum(squaresT**2)

DFM=M-1
DFE=M-N
DFT=N-1

MSM=SSM/DFM
MSE=SSE/DFE
MST=SST/DFT

R2=SSM/SST
R2_adj=1-(1-R2)*(M-1)/(M-N-1)

print("RESULT OF F-TEST(value of R2 and R2_adj respectively)")
print R2
print R2_adj
print("*******************************")
chisquared=sum(residuals**2)
Dof=M-N
chisquared_red=chisquared/Dof
p_chi2=1-scipy.stats.chi2.cdf(chisquared,Dof)
stderr_reg=scipy.sqrt(chisquared_red)
chisquare=(p_chi2,chisquared,chisquared_red,Dof,R2,R2_adj)
print("CHISQUARE TEST RESULT(p_chi2,chisquared,chisquared_red,Dof,R2,R2_adj)")
print chisquare
print("*******************************")
"""THIS IS THE SHAPIRO TEST TO ANALYZE RESIDUALS"""
w,p_shapiro=scipy.stats.shapiro(residuals)
mean_res=scipy.mean(residuals)
stddev_res=scipy.sqrt(scipy.var(residuals))
t_res=mean_res/stddev_res
p_res=1-scipy.stats.t.cdf(t_res,M-1)
print("RESULT OF SHAPIRO RESIDUAL ANALYSIS(p_res)")
print p_res
"""if p_res<0.05:
Null hypothesis is rejectedn and mean is non zero
it should be high for a good fit"""

"""THIS IS THE F-TEST ON RESIDUALS""" 
F=MSM/MSE
p_F=1-scipy.stats.f.cdf(F,DFM,DFE)
print("*******************************")
print("RESULT OF F-TEST ON RESIDUALS")
"""if p_F <0.05, Null hypothesis is rejected
THAT IS TO SAY THAT R**2>0 and atleast one of the fitting parameters >0"""

dw=stools.durbin_watson(residuals)
print("DURBIN WATSON)")
resanal=(p_shapiro,w,mean_res,p_res,F,p_F,dw)
print dw
print("*******************************")
plt.title('FREUNDLICH:PLOT OF log(x/m) vs log(C) to follow Freundlichs Isotherm')
plt.show()
