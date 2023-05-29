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
  
#initial guess surface
theta0=0.000001
thetafinal=np.pi-0.000001
dtheta=(thetafinal-theta0)/(n-1)
theta=np.linspace(theta0,thetafinal,n)

h=0.5*theta**2+2
M=1

#one black hole 

def fpsi(theta,h,M):
    return 1+M/(2*h)
def drfpsi(theta,h,M):
    return -M/(2*h**2)
def dthetafpsi(theta,h,M):
    return 0