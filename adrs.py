import math
import numpy as np
import matplotlib.pyplot as plt

#u,t = -V u,x + k u,xx  -lamda u + f

# PHYSICAL PARAMETERS
K = 0.1     #Diffusion coefficient
L = 1.0     #Domain size
Time = 20.  #Integration time

plt.figure(1)
plt.figure(2)


V=1
lamda=1

# NUMERICAL PARAMETERS
NX = 5  #Number of grid points
NT = 10000   #Number of time steps max
ifre=1000000  #plot every ifre time iterations
eps=0.001     #relative convergence ratio
niter_refinement=20      #niter different calculations with variable mesh size

errorL2=np.zeros((niter_refinement))
errorH1=np.zeros((niter_refinement))
tab_h=np.zeros((niter_refinement))

for iter in range (niter_refinement):
    NX=NX+5

    dx = L/(NX-1)                 #Grid step (space)
    dt = dx**2/(V*dx+K+dx**2)   #Grid step (time)  condition CFL de stabilite 10.4.5
    print(dx,dt)

    ### MAIN PROGRAM ###

    # Initialisation
    x = np.linspace(0.0,1.0,NX)
    T = np.zeros((NX)) #np.sin(2*np.pi*x)
    F = np.zeros((NX))
    rest = []
    RHS = np.zeros((NX))

    Tex = np.zeros((NX)) #np.sin(2*np.pi*x)
    Texx = np.zeros((NX)) #np.sin(2*np.pi*x)
    for j in range (0,NX):
        xx=j/NX
        Tex[j] = np.sin(2*math.pi*xx)*np.exp(-20*(xx-0.5)**2)
    for j in range (1,NX-1):
        Texx[j]=(Tex[j+1]-Tex[j-1])/(2*dx)  #np.cos(j*math.pi/NX)*math.pi/NX  
        Txx=(Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)  #-np.sin(j*math.pi/NX)*(math.pi/NX)**2    #
        F[j]=V*Texx[j]-K*Txx+lamda*Tex[j]
        
        
    dt = dx**2/(V*dx+2*K+abs(np.max(F))*dx**2)   #Grid step (time)  condition CFL de stabilite 10.4.5

    plt.figure(1)


    # Main loop en temps
    #for n in range(0,NT):
    n=0
    res=1
    res0=1
    while(n<NT and res/res0>eps):
        n+=1
    #discretization of the advection/diffusion/reaction/source equation
        res=0
        for j in range (1, NX-1):
            xnu=K+0.5*dx*abs(V) 
            Tx=(T[j+1]-T[j-1])/(2*dx)
            Txx=(T[j-1]-2*T[j]+T[j+1])/(dx**2)
            RHS[j] = dt*(-V*Tx+xnu*Txx-lamda*T[j]+F[j])
            res+=abs(RHS[j])

        for j in range (1, NX-1):
            T[j] += RHS[j]
            RHS[j]=0

        T[NX-1]=T[NX-2]
        #T[NX-1]=(T[NX-2]*(1-2/dx**2)+T[NX-3]/dx**2)/(1-1/dx**2)

        if (n == 1 ):
            res0=res

        rest.append(res)
    #Plot every ifre time steps
        if (n%ifre == 0 or (res/res0)<=eps):
            print(n,res)
            plotlabel = "t = %1.2f" %(n * dt)
            plt.plot(x,T, label=plotlabel,color = plt.get_cmap('copper')(float(n)/NT))
          
    print(n,res)
    plt.figure(1)
    plt.plot(x,T)
    plt.plot(x,Tex)

    # plt.xlabel(u'$x$', fontsize=26)
    # plt.ylabel(u'$T$', fontsize=26, rotation=0)
    # plt.title(u'ADRS 1D')
    # plt.legend()

    plt.figure(2)
    plt.plot(np.log10(rest/rest[0]))

    errL2=np.dot(T-Tex,T-Tex)/NX
    
    errh1=0
    h1_ex=0
    for j in range (1,NX-1):
        Txx_ex=(Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)
        h1_ex+=Txx_ex**2/NX
        errh1+=(Texx[j]-(T[j+1]-T[j-1])/(2*dx))**2/NX
       
    tab_h[iter]=1/NX
    errorL2[iter]=np.sqrt(errL2)/np.sqrt(h1_ex)
    errorH1[iter]=np.sqrt(errh1)/np.sqrt(h1_ex)
    print('norm error=',errorL2[iter],errorH1[iter])

f=plt.figure(3)
f.add_subplot(1,2, 1)
plt.plot(errorL2,label="L2")
plt.plot(errorH1,label="H1")
plt.legend()
f.add_subplot(1,2, 2)
plt.plot(np.log10(rest/rest[0]),label="Convergence in time")
plt.legend()

#%%

plt.figure(4)
# plt.plot(x,Tex, label=plotlabel,color = plt.get_cmap('copper')(float(n)/NT))

from scipy.optimize import curve_fit

#modele lineaire
def func(x, a, b):
    return a * x + b

def func_exp(x, a, b):
    return b*x**a

xdata=np.log10(tab_h)
ydata=np.log10(errorL2)
popt, pcov = curve_fit(func, xdata, ydata)

model = func(xdata,popt[0],popt[1])

print("code is order",popt[0],"constant C=",10**popt[1])
print("cov:",pcov)

#modele exponentiel
xdata=tab_h
ydata=errorL2
popt, pcov = curve_fit(func_exp, xdata, ydata)

model_exp = func_exp(xdata,popt[0],popt[1])

plt.plot(ydata,label="L2")
plt.plot(10**model,label="model")
plt.plot(model_exp,label="model_exp")
plt.legend()

print("code is order",popt[0],"constant C=",popt[1])
print("cov:",pcov)

