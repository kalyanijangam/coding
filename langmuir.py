import scipy
import matplotlib.pyplot as plt
#import statsmodels
import statsmodels.stats.stattools as stools
import scipy.optimize
import pandas as pd


filehandle=open('C:\\abs.xlsx','r')
print filehandle

'''standardization of NaOH'''
eb=0.1
ep=0.1
ewb=0.001

N_oxalic=0.2
V_oxalic=5
V_naoh=22.9

'''Normality of NaOH= (N_oxalic*V_oxalic/V_naoh)
Error(normality)/normality=((error(V_oxalic)/V_oxalic)^2-(error(V_naoh)/V_naoh)^2)^0.5 
'''

e1=V_oxalic/ep
e2=V_naoh/eb
e3=(abs((e1**2.0)-(e2**2.0)))**0.5
N_naoh=N_oxalic*(V_oxalic/V_naoh)
e_standardization=e3*N_naoh


e=pd.ExcelFile('C:\\abs.xlsx')
sheet4=e.parse(3)
V_nacl=sheet4.icol(2).real
ep=sheet4.icol(8).real
V_mixture=sheet4.icol(10).real
eb=sheet4.icol(9).real
V_Naoh=sheet4.icol(3).real
N_hcl=N_naoh*V_Naoh/V_mixture

e4=V_mixture/ep
e5=V_Naoh/eb
e6=(abs((e4**2.0)-(e5**2.0)+(e3**2.0)))**0.5

e_strength=e6*N_hcl

#gm equivalents of HCl formed=x

x=(N_hcl/1000)*V_nacl
e7=V_nacl/ep
e8=(abs((e6**2.0)-(e7**2.0)))**0.5
e_x=e8*x

#(x/m) ratio

m=sheet4.icol(1).real
q=x/m
ewb=0.001
e9=(abs((e8**2.0)-(ewb**2.0)))**0.5
e_q=e9*q
M1=(x/m)+e_q
M2=(x/m)-e_q

c1=sheet4.icol(8).real

c=c1-N_hcl
e_c=(e_strength)*c

print "Normality of standardizeedd NaOH=:\n",N_naoh
print "strength of HCl in solution=:\n",N_hcl
print "gram equivalents of HCl:\n",x
print "value of x/m is:\n",q
print "error in calculation of the strength of HCl is:\n",e_strength
print "error in calculation of x:\n",e_x
print "error in calculation of x/m:\n",e_q
print "error in calculation of c:\n",e_c
print "absolute error in the experiment is:\n", e_q

u=1/(x/m)
t=1/c

def curve(t,p):
    [A,B]=p
    u=A+(B*t)
    return u
def error(p,t,uexp):
    ucalc=curve(t,p)
    err=ucalc-uexp
    return err
def get_r2(t,u,ucalc):
    umean=scipy.average(u)
    sigmai=(((u-umean)**2)/9)**0.5
    dumean2=(u-umean)**2  
    ducalc2=(u-ucalc)**2
    r2=1-sum(ducalc2)/sum(dumean2)
    return r2,sigmai
umean=scipy.average(u)


sigmai=(((u-umean)**2)/9)**0.5

pguess=[1.050,0.587]
plsq=scipy.optimize.leastsq(error,pguess,args=(t,u))
p=plsq[0]
ucalc=curve(t,p)
j=u-ucalc

chisquare=sum((j**2)/((sigmai**2)*81))


sigmasq=sum(((u-ucalc)**2)/9)
sigmau=scipy.sqrt(sigmasq)

dumean2=(u-umean)**2  
ducalc2=(u-ucalc)**2

dsducalc2=sum(ducalc2)
dsdumean2=sum(dumean2)
r2=1-dsducalc2/dsdumean2

residual=(u-ucalc)/e_q #calculating the residuals
w,p_shapiro=scipy.stats.shapiro(residual)#shapiro wilk test
dw=stools.durbin_watson(residual)#durbin watson test

DFM=8
DFE=1


squares=(ucalc-umean)
squaresT=(u-umean)
residuals=(ucalc-u)

SSM=sum(squares**2)
SSE=sum(residuals**2)
SST=sum(squaresT**2)

MSM=SSM/DFM
MSE=SSE/DFE

F=MSM/MSE

p_F=1-scipy.stats.f.cdf(F,DFM,DFE)

print "F test value is:",p_F

print "chisquare=:\n",chisquare

print "Durbin Watson parameter: \n",dw
print "r2=:\n",r2

print "Shapiro parameter:\n",w 


print "sigmai=:\n",sigmai

print "uncertainty in data=:\n",sigmau

print "residuals are=:\n",residual

print "calculated values of log(x/m) from curve fit=:\n",ucalc

print "values of coefficients :\n", p


k=1/(x/m)

M1=ucalc+sigmau
M2=ucalc-sigmau

fig=plt.figure()
fig.show
axy=fig.add_subplot(111)
axy.title.set_text('m/x vs 1/c')
axy.plot(1/c,k,'o')
axy.plot(1/c,ucalc,label='data fit')
axy.plot(1/c,M1,label='Standard error of fit upper limit')
axy.plot(1/c,M2,label='Standard error of fit lower limit')
plt.legend()


fig=plt.figure()
fig.show
axy=fig.add_subplot(111)
axy.title.set_text('residuals vs u using fit')
axy.plot(u,residual,'*')






