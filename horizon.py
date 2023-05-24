"""
The code below was written by @author: https://github.com/DianaNtz and
solves the Apparent Horizon equation in axi symmetry for a time symmetric 
conformally flat metric. The code is an implementation of the shooting 
method which solves a boundary value problem by reducing it to an initial 
value problem.
"""
import numpy as np
import matplotlib.pyplot as plt
#Brill-Lindquist

#One black hole 
def fpsi(theta,r,M):
    return 1+M/(2*r)

def drfpsi(theta,r,M):
    return -M/(2*r**2)

def dthetafpsi(theta,r,M):
    return 0
def fr1(r1,r2,theta,M):
    return r2

def fr2(r1,r2,theta,M):
    
    if(np.abs(theta) < 10**(-8) or np.abs(theta-np.pi) < 10**(-8) ):
        factor=0
    else:
        factor=r2*(1+(r2/r1)**2)*(np.cos(theta)/np.sin(theta))
        
    value=2*r1-factor+4*r1**2*(1+(r2/r1)**2)/fpsi(theta,r1,M)*(drfpsi(theta,r1,M)
         -(1/r1**2)*r2*dthetafpsi(theta,r1,M))+3*r2**2/r1
    if(np.abs(r1) > 1.2182772912493089e+14 or np.abs(r1) < 10**(-8)):        
        return 0
    else:
        return value