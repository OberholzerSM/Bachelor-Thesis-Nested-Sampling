# -*- coding: utf-8 -*-
"""
Created on Thu May 21 09:43:31 2020

@author: samue
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import dynesty

n=100;
u=6;

wj = (2*math.pi) / 365.25;
wc = (2*math.pi) / 432.25;

X = np.zeros(n);
Y = np.zeros(n);
dX = np.zeros(n);
dY = np.zeros(n);
Time = np.zeros(n);

Xt = np.loadtxt('CODE_X.txt')
Yt = np.loadtxt('CODE_Y.txt')
dXt = np.loadtxt('CODE_dX.txt')
dYt = np.loadtxt('CODE_dY.txt')
MJD = np.loadtxt('CODE_MJD.txt')

for i in range (0,n):
    
    X[i] = Xt[9673 - n + i];
    Y[i] = Yt[9673 - n + i];
    
    if dXt[9673 - n + i] == 0:
        dX[i] = math.pow(10,-6);
    else:
        dX[i] = dXt[9673 - n + i];
        
    if dYt[9673 - n + i] == 0:
        dY[i] = math.pow(10,-6);
    else:
        dY[i] = dYt[9673 - n + i];
    
    Time[i] = MJD[9673 - n + i];
    
def ptform(uniform):
    
    z = np.zeros(u)
    for i in range(0,u):
        z[i] = 2*uniform[i] - 1
    
    return z;


summe = 0;
for i in range(0,n):
    
    summe += -1* ( math.log(dX[i]) + math.log(dY[i]) + (math.pow(X[i],2) / (2*dX[i])) + (math.pow(Y[i],2) / (2*dY[i]))  );

summe += -n*math.log(2*math.pi)

def loglike(theta):
    
    def x(t):
        return theta[0]*math.cos(wj*t) - theta[1]*math.sin(wj*t) + theta[2]*math.cos(wc*t) - theta[3]*math.sin(wc*t) + theta[4];
        
    def y(t):
        return -theta[0]*math.cos(wj*t) - theta[1]*math.sin(wj*t) - theta[2]*math.cos(wc*t) - theta[3]*math.sin(wc*t) + theta[5];

    lnL = summe;
    for i in range(0,n):
        lnL += ( -math.pow( x(Time[i]) , 2) / (2*dX[i]) ) + ( (x(Time[i]) * X[i]) / dX[i] ) + ( -math.pow( y(Time[i]) , 2) / (2*dY[i]) ) + ( (y(Time[i]) * Y[i]) / dY[i] ) ;
    
    
    return lnL
    
sampler = dynesty.NestedSampler(loglike, ptform, u);
sampler.run_nested();
results = sampler.results;

samples = results.samples
lnZ = results.logz[len(samples)-1];
lnw = results.logwt;

p = np.zeros(len(samples));
P = np.zeros(len(samples));
Rxj = np.zeros(len(samples));
Ryj = np.zeros(len(samples));
Rxc = np.zeros(len(samples));
Ryc = np.zeros(len(samples));
x0 = np.zeros(len(samples));
y0 = np.zeros(len(samples));


for i in range(0,len(samples)):
    
    p[i] = math.exp( lnw[i] - lnZ  );
    P[i] = 100*math.exp( lnw[i] - lnZ  );
    Rxj[i] = samples[i,0];
    Ryj[i] = samples[i,1];
    Rxc[i] = samples[i,2];
    Ryc[i] = samples[i,3];
    x0[i] = samples[i,4];
    y0[i] = samples[i,5];
    
    

plt.figure(1)
plt.plot(Rxj,P)
plt.figure(2)
plt.plot(Ryj,P)
plt.figure(3)
plt.plot(Rxc,P)
plt.figure(4)
plt.plot(Ryc,P)
plt.figure(5)
plt.plot(x0,P)
plt.figure(6)
plt.plot(y0,P)

plt.figure(7)
plt.plot(Rxj)
plt.figure(8)
plt.plot(Ryj)
plt.figure(9)
plt.plot(Rxc)
plt.figure(10)
plt.plot(Ryc)
plt.figure(11)
plt.plot(x0)
plt.figure(12)
plt.plot(y0)

plt.show();

print("Mittelwerte:");
print( sum(p*Rxj) , sum(p*Ryj), sum(p*Rxc), sum(p*Ryc) ,sum(p*x0) ,sum(p*y0));
print(" ");
print("Zweite Momente:");
print( sum(p*Rxj**2) , sum(p*Ryj**2), sum(p*Rxc**2), sum(p*Ryc**2) ,sum(p*x0**2) ,sum(p*y0**2));