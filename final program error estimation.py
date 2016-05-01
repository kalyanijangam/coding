# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 21:16:36 2016

@author: USER
"""

import scipy,numpy
import scipy.optimize, scipy.stats
import numpy.random
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels
import statsmodels.stats
import statsmodels.stats.stattools as stools

#import statsmodels.formula.api as sm


plt.style.use("ggplot")

'''formatting of axes of graphs'''
def formataxis(ax):
    ax.xaxis.label.set_fontname('georgia')
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontname('georgia')
    ax.yaxis.label.set_fontsize(12)
    ax.title.set_fontname('georgia')
    ax.title.set_fontsize(12)
    
    
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)

'''calculation of standard error of fit'''
# M= number of data points and N= no of parameters present in function
def get_stderr_fit(f,xdata,popt,pcov,dict_data):
    Y= f(xdata,popt,dict_data)
    listdY=[]
    for i in xrange(len(popt)):
        p=popt[i]
        dp=abs(p)/1e6+1e-20
        popt[i]+=dp
        Yi= f(xdata,popt,dict_data)
        dY=(Yi-Y)/dp
        listdY.append(dY)
        popt[i]-= dp
    listdY= scipy.array(listdY)
    left= scipy.dot(listdY.T,pcov)
    right= scipy.dot(left,listdY)
    sigma2y= right.diagonal()
    mean_sigma2y=scipy.mean(right.diagonal())
    M=Xdata.shape[1]
    N=len(popt)
    avg_stddev_data= scipy.sqrt(M*mean_sigma2y/N)
    sigmay=scipy.sqrt(sigma2y)
    #sigmay : the standard error of fit and is a function of x
    return sigmay,avg_stddev_data
    popt=numpy.array([0.0]*11)
#popt=optimized parameters


def fitdata(f,Xdata,Ydata,Errdata,pguess,dict_data,ax=False,ax2=False):
 
     def error(p, Xdata, Ydata, Errdata, dict_data):
        Y=f(Xdata, p,dict_data)
        residuals= (Y-Ydata)/Errdata
        
        return residuals
     res = scipy.optimize.leastsq(error, pguess, args=(Xdata, Ydata, Errdata, dict_data), full_output=1)
     (popt, pcov, infodict, errmsg, ier) = res
     perr=scipy.sqrt(scipy.diag(pcov))
     M= len(Ydata)
     N=len(popt)
        #residuals
     Y=f(Xdata,popt,dict_data)
     residuals=(Y-Ydata)/Errdata
     meanY=scipy.mean(Ydata)
     squares=(Y-meanY)/Errdata
     squaresT=(Ydata-meanY)/Errdata
    
     SSM=sum(squares**2) #Corrected Sum of Squares
     SSE=sum(residuals**2) #Sum of Squares of Errors
     SST=sum(squaresT**2) #Total corrected sum of squares
    
     DFM=N-1 #Degrees of freedom for model
     DFE=M-N #Degrees of freedom for error
    # DFT=M-1 #Degrees of freedom total
    
     MSM=SSM/DFM #Mean squares for model (explained variance)
     MSE=SSE/DFE #Mean squares for Error (should be small wrt MSM) Unexplained Variance
     #MST=SST/DFT #Mean squares for total
    
     R2=SSM/SST #proportion of explained variance
     R2_adj=1-(1-R2)*(M-1)/(M-N-1) #Adjusted R2
        
    #t-test to see if parameters are different from zero
     t_stat=popt/perr #t-statistic for popt different from zero'
     t_stat=t_stat.real
     p_p=1.0-scipy.stats.t.cdf(t_stat,DFE) #should be low for good fit.
     z=scipy.stats.t(M-N).ppf(0.95)
     p95=perr*z
    #Chisquared Analysis on Residuals
     chisquared=sum(residuals**2)
     degfreedom=M-N
     chisquared_red=chisquared/degfreedom
     p_chi2=1.0-scipy.stats.chi2.cdf(chisquared,degfreedom)
     stderr_reg=scipy.sqrt(chisquared_red)
     chisquare=(p_chi2,chisquared, chisquared_red, degfreedom,R2,R2_adj)
    
    #Analysis of Residuals
     w,p_shapiro=scipy.stats.shapiro(residuals)
     mean_res=scipy.mean(residuals)
     stddev_res=scipy.sqrt(scipy.var(residuals))
     t_res=mean_res/stddev_res #t-statistic to test that mean_res is zero.
     p_res=1.0-scipy.stats.t.cdf(t_res,M-1)
        #if p_res <0.05, null hypothesis rejected and mean is non-zero.
        #Should be high for good fit.
    #F-test on residuals
     F=MSM/MSE #explained variance/unexplained . Should be large
     p_F=1.0-scipy.stats.f.cdf(F,DFM,DFE)
        #if p_F <0.05n, null-hypothesis is rejected
        #i.e. R^2>0 and at least one of the fitting parameters >0.
     dw=stools.durbin_watson(residuals)
     resanal=(p_shapiro,w,mean_res,F,p_F,dw)
     print "absolute error: \n",residuals
     
     if ax:
         formataxis(ax)
         ax.plot(Ydata,Y,'ro')
         ax.errorbar(Ydata,Y,yerr=Errdata,fmt='.')
         Ymin,Ymax=min((min(Y),min(Ydata))),max((max(Y),max(Ydata)))
         ax.plot([Ymin,Ymax],[Ymin,Ymax],'b')
        
         ax.xaxis.label.set_text('Data')
         ax.yaxis.label.set_text('Fitted')
        
         sigmay,avg_stddev_data=get_stderr_fit(f,Xdata,popt,pcov,dict_data)
         Yplus=Y+sigmay
         Yminus=Y-sigmay
         print "confidence range: \n"
         print "Y + standard deviation:  \n", Yplus
         print "--------------------------------------------------------"
         print "Y: \n", Y
         print "--------------------------------------------------------"
         print "Y - standard deviation : \n", Yminus
         ax.plot(Y,Yplus,'k',alpha=1.0,linewidth=0.8) #,linestyle='--'
         ax.plot(Y,Yminus,'k',alpha=1.0,linewidth=0.8)
         ax.fill_between(Y,Yminus,Yplus,facecolor='cyan',alpha=1.0)
         titletext='Parity Plot for Fitted Data \n'
         titletext+=r'$r^2$=%5.2f, $r^2_{adjusted}$=%5.2f, '
         titletext +=r' $\sigma_{exp}$=%5.2f,$\chi^2_{\nu}$=%5.2f, $p_{\chi^2}$=%5.2f, '
         titletext +=r'$\sigma_{err}^{reg}$=%5.2f'
         ax.title.set_text(titletext%(R2,R2_adj, avg_stddev_data, chisquared_red,  p_chi2, stderr_reg))
         ax.figure.canvas.draw()
         
        
     if ax2:#Test for homoscedasticity
        formataxis(ax2)
        ax2.plot(Y,residuals,'ro')
        
        ax2.xaxis.label.set_text('Fitted Data')
        ax2.yaxis.label.set_text('Residuals')
        
        titletext='Analysis of Residuals\n'
        titletext+=r'mean=%5.2f, $p_{res}$=%5.2f, $p_{shapiro}$=%5.2f, $Durbin-Watson$=%2.1f'
        titletext+='\nF=%5.2f, $p_F$=%3.2e'
        ax2.title.set_text(titletext%(mean_res, p_res, p_shapiro, dw , F, p_F))
        
        ax2.figure.canvas.draw()
        
     return popt,pcov,perr,p95,p_p,chisquare,resanal
         
def import_data(xlfile,sheetname):
    df=pd.read_excel(xlfile,sheetname=sheetname)
    return df
    
def prepare_data(df,Criterion,Predictors,Error=False):
    Y=scipy.array(df[Criterion])
    if Error:
        Errdata=scipy.array(df[Error])
    else:
        Errdata=scipy.ones(len(Y))
    Xdata=[]
    for X in Predictors:
        X=list(df[X])
        Xdata.append(X)
    Xdata=scipy.array(Xdata)
    
    return Xdata,Y,Errdata
    
if  __name__ =="__main__":
    fig=plt.figure()
    ax=fig.add_subplot(111)
    fig.show()
    
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    fig2.show()
    
# Make arbitrary function of three variables
    def f(X,p,dict_data):
#        a=dict_data['a']
#        b=dict_data['b']
        (x,x)=X
        Y=p[0]+p[1]*x**2+p[2]*x
        return Y
   
    
    
    
    #Get data from excel file using pandas
    dict_data=[]
    df=import_data('teabag infusion data.xlsx','Data')
    Xdata,Ydata,Errdata = prepare_data(df,'y',('x','x'),Error= 'err')
    #Initial Guess
    N=3
    pguess=N*[0.0]
    

    
    popt,pcov,perr,p95,p_p,chisquare,resanal=fitdata(f, Xdata, Ydata, Errdata, pguess,dict_data, ax=ax, ax2=ax2)
    print "--------------------------------------------------------"
    print "optimum parameters of function: ",popt
    print "the standard-deviation of the optimized parameters: \n", perr
    print "covariance matrix of the fit :\n", pcov
#    from tabulate import tabulate
#    table = [[popt],[perr]]
#    print tabulate(table,headers=["optimum parameters", "standard deviation"], tablefmt="fancy_grid")
    
    
        
        
        
        