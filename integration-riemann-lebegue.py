#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:36:14 2024

@author: bijan.mohammadi@umontpellier.fr
"""
from scipy.integrate import quad
import math
import numpy as np
import matplotlib.pyplot as plt

def fct(x, a, b, c):
    return a*x**2 + b*x + c *np.sin(4*np.pi*x) + 10*np.exp(-100*(x-0.5)**2)

def fct_xx(x, a, b, c):
    return 2*a - c*(4*np.pi)**2*np.sin(4*np.pi*x) + 10*40000*(x-0.5)*np.exp(-100*(x-0.5)**2)

#reference
L = 1
a = 2
b = 1
c=3
I = quad(fct, 0, L, args=(a,b,c))
print(I)

#%%
npt=1000
x=0
hh=L/(npt-1)
for i in range(npt):
    x=(i-1)*hh
    plt.plot(x,fct(x,a,b,c),".")
plt.xlabel("x")
plt.xlabel("F(x)")
plt.legend()
plt.show()

#%%
#integrer en Riemann avec le nombre de points d'integration croissant
#h constant

h=np.zeros(100)
IR=np.zeros(100)
for npt in range(3,100):
    h[npt]=L/(npt-1)
    for j in range(npt):
        x=(j-1)*h[npt]
        IR[npt]+=h[npt]*fct(x,a,b,c)
        
#%%
#Integrer en Lebesgues: h donne par formule (2) du poly

hmin=0.01
hmax=0.5
epsilon=1

itermax_lebesgue=100
nptL=np.zeros(itermax_lebesgue).astype(int)
IL=np.zeros(itermax_lebesgue)
for npt in range(itermax_lebesgue):
    epsilon*=0.9
    x=0
    u=fct(x, a, b, c)
    while(x<L):
        uxx=fct_xx(x, a, b, c)
        metric=min(max(abs(uxx)/epsilon,1/hmax**2),1/hmin**2)
        hloc=min(np.sqrt(1./metric),L-x)
        #print(hloc**2*metric)
        u0=u
        x+=hloc
        u=fct(x, a, b, c)
        IL[npt]+=hloc*(u+u0)/2
        nptL[npt]+=1
    print(epsilon,nptL[npt])

plt.plot(np.log10(abs(IR[3:]-I[0])),label="Riemann")        #1/h[3:],
plt.plot(np.log10(abs(IL-I[0])),"-",label="Lebesgues")
plt.xlabel("ITER REFINEMENT")
plt.legend()
plt.show()

plt.plot(nptL,label="NB PT BY LEBESGUE")
plt.xlabel("ITER LEBESGUE")
plt.legend()
plt.show()

