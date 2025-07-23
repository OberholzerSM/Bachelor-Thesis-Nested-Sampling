# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 14:11:57 2020

@author: samue
"""

import math
import numpy as np
import dynesty
from dynesty import plotting as dyplot

n=21;
u=2;

X = np.zeros(n);
Y = np.zeros(n);

for i in range(0,n):
    
    a = i - 10;
    
    X[i] = a;
    Y[i] = np.random.normal((2*a)-1,0.5);
    
def ptform(u):
    return [2*u[0]+1,2*u[1]-2];

def loglike(theta):
    
    def y(x):
        return theta[0]*x + theta[1];

    lnL = 0;
    for i in range(0,n):
        lnL += math.pow(y(X[i])-Y[i],2);
    
    lnL = -0.5 * lnL - ((n/2)*math.log(0.5*math.pi));
    return lnL
    
sampler = dynesty.NestedSampler(loglike, ptform, u);
sampler.run_nested();
results = sampler.results;

rfig, raxes = dyplot.runplot(results)

fig, axes = dyplot.traceplot(results, truths=np.zeros(u),
                             truth_color='black', show_titles=True,
                             trace_cmap='viridis', connect=True,
                             connect_highlight=range(5))

fg, ax = dyplot.cornerplot(results, color='blue', truths=np.zeros(u),
                           truth_color='black', show_titles=True)