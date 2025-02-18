# -*- coding: utf-8 -*-
"""
ref:
LBM mesoscale modelling of porous media
A. C. M. Sousa1,2 & A. Nabovati1
1Department of Mechanical Engineering, University of New Brunswick,
Fredericton, NB, Canada
2Department of Mechanical Engineering, University of Aveiro, Aveiro,
Portugal

ref:
Slow flow through a periodic array of spheres
Author links open overlay panelA.S.SanganiA.Acrivos

ref:
LBflow: An extensible lattice Boltzmann framework for the simulation of geophysical flows. Part II: usage and validation☆
Author links open overlay panelE.W.Llewellin

Note: this semi-analytical solution works best we porosity is in range of 0.3-0.5
"""
import numpy as np
r = 35e-9
L = 80e-9
cmax= 3**0.5 * np.pi/8
c = 8 * np.pi * r**3/(3*L**3)
x = (c/cmax)**(1/3)
dst = 0 
alpha = [0.1000000e1,
         0.1575834e01,
         0.2483254e01,
         0.3233022e01,
         0.4022864e01,
         0.4650320e01,
         0.5281412e01,
         0.5826374e01,
         0.6258376e01,
         0.6544504e01,
         0.6878396e01,
         0.7190839e01,
         0.7268068e01,
         0.7304025e01,
         0.7301217e01,
         0.7236410e01,
         0.7298014e01,
         0.7369849e01,
         0.7109497e01,
         0.6228418e01,
         0.5235796e01,
         0.4476874e01,
         0.3541982e01,
         0.2939353e01,
         0.3935484e01,
         0.5179097e01,
         0.3959872e01,
         0.2227627e01,
         0.3393390e01,
         0.4491369e01,
         0.2200686e01]
for i in range(31):
    dst+= alpha[i] * x**i
nsphere =  8 * 1/8 + 1
k =L**3/(6 * nsphere *np.pi * r *dst)

vsphere = 4/3 * np.pi* r**3
nsphere =  8 * 1/8 + 1
phi =  1- ((nsphere * vsphere)/(L**3))
print(phi,1-phi, k)
phe = np.log(phi/(1-phi))
k = r**2 * np.exp(-4.6254 + 2.2648*phe - 0.2093*phe**2 + 0.118*phe**3)
print(k)
#print( 2 *r**2/(9*phi*dst))