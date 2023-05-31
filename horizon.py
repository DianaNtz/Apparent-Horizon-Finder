"""
The code below was written by @author: Diana Nitzschke (https://github.com/DianaNtz) and
solves the Apparent Horizon equation in axi symmetry for a time symmetric 
conformally flat metric. The code is an implementation of the flow 
method which solves an elliptic PDE by transforming it into a hyperbolic equation. 
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
filenames = []

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

#Brill-Lindquist
M=1
s=0.75
#one black hole 
"""
def fpsi(theta,h,M):
    return 1+M/(2*h)
def drfpsi(theta,h,M):
    return -M/(2*h**2)
def dthetafpsi(theta,h,M):
    return 0
"""
#two black holes with same mass
"""
def fpsi(theta,r,M):
    return 1+0.5*M/np.sqrt(r**2+2*s*r*np.cos(theta)+(s)**2)+0.5*M/np.sqrt(r**2-2*s*r*np.cos(theta)+(s)**2)

def drfpsi(theta,r,M):
    return -0.5*M/np.sqrt(r**2+2*s*r*np.cos(theta)+(s)**2)**3*(r+s*np.cos(theta))-0.5*M/np.sqrt(r**2-2*s*r*np.cos(theta)+(s)**2)**3*(r-s*np.cos(theta))

def dthetafpsi(theta,r,M):
    return 0.5*M/np.sqrt(r**2+2*s*r*np.cos(theta)+(s)**2)**3*(s*r*np.sin(theta))-0.5*M/np.sqrt(r**2-2*s*r*np.cos(theta)+(s)**2)**3*(s*r*np.sin(theta))
"""
#two black holes with different mass

def fpsi(theta,r,M):
    return 1+0.5*M*0.8/np.sqrt(r**2+2*0*r*np.cos(theta)+(0)**2)+0.5*M*0.2/np.sqrt(r**2-2*0.65*r*np.cos(theta)+(0.65)**2)

def drfpsi(theta,r,M):
    return -0.5*M*0.8/np.sqrt(r**2+2*0*r*np.cos(theta)+(0)**2)**3*(r+0*np.cos(theta))-0.5*M*0.2/np.sqrt(r**2-2*0.65*r*np.cos(theta)+(0.65)**2)**3*(r-0.65*np.cos(theta))

def dthetafpsi(theta,r,M):
    return 0.5*M*0.8/np.sqrt(r**2+2*0*r*np.cos(theta)+(0)**2)**3*(0*r*np.sin(theta))+0.5*M*0.2/np.sqrt(r**2-2*0.65*r*np.cos(theta)+(0.65)**2)**3*(-0.65*r*np.sin(theta))

#expansion THETA
def THETA(h,theta):
    d1=d1theta(h,theta)
    d2=d2theta(h,theta)
    Fpsi=fpsi(theta,h,M)
    C=1/np.sqrt(1+d1**2/h**2)
    T=C**3/(Fpsi**2*h**2)*(-d2+2*h-d1*(1+(d1/h)**2)*(np.cos(theta)/np.sin(theta))+4*h**2*(1+(d1/h)**2)/Fpsi*(drfpsi(theta,h,M)-(1/h**2)*d1*dthetafpsi(theta,h,M))+3*d1**2/h)
    return -T*C*Fpsi**(-2)

#integration method
dlambda=0.025*dtheta
for k in range(0,60000*6+1):
    k1=dlambda*THETA(h,theta)
    h=h+k1
    if(k%1500==0):
        print(k)
        
        phi = np.linspace(0, 2.0*np.pi, 200)
        T,P=np.meshgrid(theta, phi)
        
        fig = plt.figure(constrained_layout=True,figsize=(6,6))
        ax = fig.add_subplot(projection='3d')
        ax.plot_surface(h*np.sin(T)*np.cos(P),h*np.sin(T)*np.sin(P),h*np.cos(T),  rstride=1, cstride=1,cmap='seismic')
        ax.view_init(30,90)
        #ax.set_xlabel("x",fontsize= 19)
        #ax.set_ylabel("y",fontsize= 19)
        #ax.set_zlabel("z",fontsize= 19)
        ax.zaxis.set_tick_params(labelsize=18,pad=7)
        ax.yaxis.set_tick_params(labelsize=18,pad=7)
        ax.xaxis.set_tick_params(labelsize=18)
        stringtitle="λ=".__add__(str(round(k,0))).__add__("dλ")
        plt.title(stringtitle,fontsize=20,x=0.5, y=0.95)
        filename ='bla{0:.0f}.png'.format(int(k/1500))
        #append file name to the list filename
        filenames.append(filename)    
        #save the plot
        plt.savefig(filename,dpi=100)
        plt.close()
        #plt.show()
   
        
#build the gif
with imageio.get_writer('figures/twoblackholesdifferentmassneu.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)       
#remove saved figures 
for filename in set(filenames):
    os.remove(filename)  
   
phi = np.linspace(0, 2.0*np.pi, 200)
T,P=np.meshgrid(theta, phi)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')
ax.plot_surface(h*np.sin(T)*np.cos(P),h*np.sin(T)*np.sin(P),h*np.cos(T),  rstride=1, cstride=1,cmap='seismic')
ax.view_init(30,90)
plt.show()
    