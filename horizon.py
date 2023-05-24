"""The code below was written by @author: https://github.com/DianaNtz and
solves the Apparent Horizon equation in axi symmetry for a time symmetric 
conformally flat metric. The code is an implementation of the shooting 
method which solves a boundary value problem by reducing it to an initial 
value problem."""
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
