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
#integration method
def shooting(r0,M):
    ntheta=600
    
    r10=r0
    r20=0
    
    theta0=0
    thetafinal=np.pi
    dtheta=(thetafinal-theta0)/(ntheta-1)
    
    
    theta=np.zeros(ntheta)
    r1=np.zeros(ntheta)
    r2=np.zeros(ntheta)
    
    thetan=theta0
    r1n=r10
    r2n=r20
    for j in range(0,ntheta):
        theta[j]=thetan
        r1[j]=r1n
        r2[j]=r2n
        
        k1r1=dtheta*fr1(r1n,r2n,thetan,M)
        k1r2=dtheta*fr2(r1n,r2n,thetan,M)
        
        k2r1=dtheta*fr1(r1n+0.5*k1r1,r2n+0.5*k1r2,thetan+0.5*dtheta,M)
        k2r2=dtheta*fr2(r1n+0.5*k1r1,r2n+0.5*k1r2,thetan+0.5*dtheta,M)
        
        r1n=r1n+k2r1
        r2n=r2n+k2r2
        
        
        thetan=thetan+dtheta
    return theta,r1,r2[-1]
#bisection method
def find(ra1,ra2,M):
    i=0
    tol=0.00001
    N=100
    while(i<N):
         theta1,r1,value1=shooting(ra1,M)
         if(np.abs(value1)<tol):
             print(ra1)
             return theta1,r1
         theta2,r2,value2=shooting(ra2,M)
         if(np.abs(value2)<tol):
             print(ra2)
             return theta2,r2
         
         if(value1*value2 < tol):
             ra3=(ra1+ra2)/2
             theta3,r3,value3=shooting(ra3,M)
            
             if(np.abs(value3)<tol):
                print(ra3)
                return theta3,r3
             if(value1*value3 < 0):
                ra2=ra3
             if(value2*value3 < 0):
                ra1=ra3
         else:
             print("try other initial values for your search!") 
             return theta1*0,r1*0
            
         i=i+1
    print("Reached maximum number of iterations N={0:.0f}!".format(N))
    theta1,r1,value1=shooting(ra1,M)
    return theta1*0,r1*0