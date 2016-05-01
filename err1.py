#Experiment is single stage adsorption
#plotting values of log(x/m) vs logC
import numpy as np
import matplotlib.pyplot as plt
import scipy
data = np.genfromtxt('ads.txt', delimiter=',')#importing data from notepad

x = data[:,][:,0]
y =data[:,][:,1]
a= np.vstack([x, np.ones(len(x))]).T#using the least square method we write the data in form y=a.p where p=[m,c]
m, c=np.linalg.lstsq(a,y)[0]
print "Slope = ",m   #slope and intercept of graph of log (x/m) vs log c
print "Intercept = ",c

n=1/m;
k=scipy.exp(c);

print "value of K = ",k
print "value of n = ",n

plt.plot(x, y, 'o', label='Original data', markersize=10)
plt.plot(x,m*x+c,'r',label='Fitted line', markersize=10)
plt.xlabel('log C')
plt.ylabel('log (x/m)')
plt.legend()
plt.show()