"""
The code below was written by @author: https://github.com/DianaNtz and
solves the Apparent Horizon equation in axi symmetry for a time symmetric 
conformally flat metric. The code is an implementation of the flow 
method which solves an elliptic PDE by transforming it into a hyperbolic equation. 
"""
import numpy as np
import matplotlib.pyplot as plt

n=100

def d2theta(D,theta):
    dtheta=theta[1]-theta[0]
    Dxx=np.zeros((n), dtype='double')
    for j in range(0,n):
                 if(j!=0 and j!=n-1):
                     Dxx[j]=(D[j+1]-2*D[j]+D[j-1])/(dtheta**2)
    Dxx[0]=(D[1]-2*D[0]+D[n-1])/(dtheta**2)
    Dxx[n-1]=(D[0]-2*D[n-1]+D[n-2])/(dtheta**2)
    return Dxx

def d1theta(D,theta):
    dtheta=theta[1]-theta[0]
    Dx=np.zeros((n), dtype='double')
    for j in range(0,n):
                 if(j!=0 and j!=n-1):
                     Dx[j]=(D[j+1]-D[j-1])/(2*dtheta)
    return Dx  

