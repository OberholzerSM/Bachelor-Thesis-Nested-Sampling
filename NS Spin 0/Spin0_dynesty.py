# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 12:14:24 2020

@author: samue
"""

import math
import numpy as np
import dynesty

N = 10
d = 3*N

V = 0.25
T = 100.0

v_max = 3500.0

m = 4.002602*1.66053906661 * math.pow(10.0,-27);
kB = 1.380649 * math.pow(10.0,-23);
hbar = (6.62607015 * math.pow(10.0,-34))/(2 * math.pi);

a = ( math.pow(math.pi,2) * math.pow(hbar,2) ) / ( 2.0 * m * kB * T * math.pow(V,1.0/3.0) )

k_max = ( math.pow(V,1.0/3.0) * v_max ) / (math.pi * hbar)

def prior(u):
    
    x = np.zeros(d)
    
    for j in range(0,d):
        
        x[j] = float( 1 + math.floor( k_max * u[j] ) )
        
    return x

def lnL(theta):
    
    k2 = 0
    
    for j in range(0,d):
        
        k2 = k2 + math.pow( theta[j] ,2);
    
    return float(d) * math.log(k_max) - math.log( float( math.factorial(d/3) ) ) - a * k2

sampler = dynesty.NestedSampler(lnL, prior, d)
sampler.run_nested()
sresults = sampler.results

lnZ = np.max( sresults.logz )