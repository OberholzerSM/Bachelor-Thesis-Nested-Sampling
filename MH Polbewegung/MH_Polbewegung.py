# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:15:20 2020

@author: samue
"""
import math as m
import numpy as np
import dynesty as dyn
import matplotlib.pyplot as plt

#Konstanten
wj = (2*m.pi) / 365.25;
wc = (2*m.pi) / 432.25;


#Daten
X = np.loadtxt("CODE_X.txt")
Y = np.loadtxt("CODE_Y.txt")
MJD = np.loadtxt("CODE_MJD.txt")
dX = np.loadtxt("CODE_dX.txt")
dY = np.loadtxt("CODE_dY.txt")

for i in range(0,9673):
    if dX[i]==0:
        dX[i] = m.pow(10,-6)
        
    if dY[i]==0:
        dY[i] = m.pow(10,-6)

#Nutze nur die letzten n Daten
n = 1000;

temp = np.zeros(n);

for i in range(0,n):
    temp[i] = MJD[i + (9673-n)];
MJD = temp;

for i in range(0,n):
    temp[i] = X[i + (9673-n)];
X = temp;

for i in range(0,n):
    temp[i] = Y[i + (9673-n)];
Y = temp;

for i in range(0,n):
    temp[i] = dX[i + (9673-n)];
dX = temp;

for i in range(0,n):
    temp[i] = dY[i + (9673-n)];
dY = temp;

#Prior und Likelihood
def prior_transform(u):

    return 0.2*u

def lnL(theta):
    
    def x(t):
        res = theta[0]*m.cos(wj*t) - theta[1]*m.sin(wj*t) + theta[2]*m.cos(wc*t) - theta[3]*m.sin(wc*t) + theta[4];
        return res

    def y(t):
        res = -theta[0]*m.sin(wj*t) - theta[1]*m.cos(wj*t) - theta[2]*m.sin(wc*t) - theta[3]*m.cos(wc*t) + theta[5];
        return res
    
    res = 0;
    for i in range(0,n):
        res += (-1)*m.log( 2*(m.pi)*m.pow(dX[i],2)*m.pow(dY[i],2) ) - 0.5*m.pow(X[i]-x(MJD[i]),2)*m.pow(dX[i],-2) - 0.5*m.pow(Y[i]-y(MJD[i]),2)*m.pow(dY[i],-2);
        
    return res

#dynesty
sampler = dyn.NestedSampler(lnL, prior_transform, 6);
sampler.run_nested();
sresults = sampler.results.samples;

sresults = np.transpose(sresults);

plt.subplot(311); 
for i in range(0,6):
    plt.plot(sresults[i]);

nsample = round(sresults.size/6);  
theta = np.zeros(6);
for i in range(0,nsample):
    for j in range(0,6):
        theta[j] += sresults[j,i];

theta = theta/nsample;
print(theta);

x = np.zeros(n);
y = np.zeros(n);
RMS = 0;
for i in range(0,n):
    x[i] = theta[0]*m.cos(wj*MJD[i]) - theta[1]*m.sin(wj*MJD[i]) + theta[2]*m.cos(wc*MJD[i]) - theta[3]*m.sin(wc*MJD[i]) + theta[4];
    y[i] = -theta[0]*m.sin(wj*MJD[i]) - theta[1]*m.cos(wj*MJD[i]) - theta[2]*m.sin(wc*MJD[i]) - theta[3]*m.cos(wc*MJD[i]) + theta[5];
    
    RMS += m.pow(X[i]-x[i],2) + m.pow(Y[i]-y[i],2);

plt.subplot(312); 
plt.plot(MJD,x);
plt.plot(MJD,X);

plt.subplot(313); 
plt.plot(MJD,y);
plt.plot(MJD,Y);

plt.show;

#RMS
RMS = RMS / (2*n - 6);
RMS = m.sqrt(RMS);
print(RMS); 
