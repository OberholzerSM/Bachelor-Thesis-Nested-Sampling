# -*- coding: utf-8 -*-
"""
Created on Sat May 16 11:10:38 2020

@author: samue
"""

import math
import numpy as np
import dynesty
from dynesty import plotting as dyplot

n = 1000;
k = 502;
u = 1;
    
def ptform(u):
    return [u[0]];

def loglike(theta):
    
    lnL = math.log( math.factorial(n) / (math.factorial(k) * math.factorial(n-k)) ) + k*math.log(theta[0]) + (n-k)*math.log(1-theta[0]);
    
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