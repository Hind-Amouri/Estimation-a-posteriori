import math
import numpy as np
import matplotlib.pyplot as plt


def fex(NX,dx,time):
    F = np.zeros((NX))
    Tex = np.zeros((NX)) #np.sin(2*np.pi*x)
    Text = np.zeros((NX)) #np.sin(2*np.pi*x)
    Texx = np.zeros((NX)) #np.sin(2*np.pi*x)
    for j in range (1,NX-1):
        v=(np.exp(-1000*((j-NX/3)/NX)**2)+np.exp(-10*np.exp(-1000*((j-NX/3)/NX)**2)))\
            *np.sin(5*j*math.pi/NX)
        Tex[j] = np.sin(4*math.pi*time)*v
        Text[j] = 4*math.pi*np.cos(4*math.pi*time)*v
        
    for j in range (1,NX-1):
        Texx[j]=(Tex[j+1]-Tex[j-1])/(2*dx)  #np.cos(j*math.pi/NX)*math.pi/NX  
        Txx=(Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)  #-np.sin(j*math.pi/NX)*(math.pi/NX)**2    #
        F[j]=V*Texx[j]-K*Txx+lamda*Tex[j]+Text[j]
    return F,Tex,Texx

#u,t = -V u,x + k u,xx  -lamda u + f

# PHYSICAL PARAMETERS
K = 0.1     #Diffusion coefficient
L = 1.0     #Domain size
Time = 1.  #Integration time


V=1
lamda=1

# NUMERICAL PARAMETERS
NX = 5  #Number of grid points
NT = 10000   #Number of time steps max
ifre=100  #plot every ifre time iterations
eps=0.001     #relative convergence ratio
niter_refinement=20      #niter different calculations with variable mesh size

irk_max=4
alpha=np.zeros(irk_max)
for irk in range(irk_max):
    alpha[irk]=1/(irk_max-irk)
    #print(alpha[irk])
# if(irk_max==3):
#     alpha[0]=0.333
#     alpha[1]=0.5
#     alpha[2]=1

error=np.zeros((niter_refinement))

NX_tab=[]
Err_tab1=[]
Err_tab2=[]

for iter in range (niter_refinement):
    NX=NX+3
    NX_tab.append(NX)
    
    dx = L/(NX-1)                 #Grid step (space)
    dt = dx**2/(V*dx+K+dx**2)   #Grid step (time)  condition CFL de stabilite 10.4.5
    print("Nbre points in space, Time step:",dx,dt)

    ### MAIN PROGRAM ###

    # Initialisation
    x = np.linspace(0.0,1.0,NX)
    T = np.zeros((NX)) #np.sin(2*np.pi*x)
    F = np.zeros((NX))
    rest = []

    plt.figure(1)


    # Main loop en temps
    #for n in range(0,NT):
    n=0
    res=1
    res0=1
    time=0
    time_total=1
    time_tab=[]
    while(time<time_total): #n<NT and res/res0>eps):
        n+=1
        F,Tex,Texx=fex(NX,dx,time)
                        
        dt = dx**2/(V*dx+2*K+abs(np.max(F))*dx**2)   #Grid step (time)  condition CFL de stabilite 10.4.5
        time+=dt
        time_tab.append(time)
        
        T0=T.copy()

        for irk in range(irk_max):
        #discretization of the advection/diffusion/reaction/source equation
            res=0
            for j in range (1, NX-1):
                xnu=K+0.5*dx*abs(V) 
                Tx=(T[j+1]-T[j-1])/(2*dx)
                Txx=(T[j-1]-2*T[j]+T[j+1])/(dx**2)
                RHS = dt*(-V*Tx+xnu*Txx-lamda*T[j]+F[j])
                res+=abs(RHS)
                T[j] = T0[j] + RHS*alpha[irk]

        if (n == 1 ):
            res0=res
        rest.append(res)
    #Plot every ifre time steps
        if (n%ifre == 0 or (res/(res0+1.e-10))<eps):
            print("iteration, residual:",n,res)
            plotlabel = "t = %1.2f" %(n * dt)
            plt.plot(x,T, label=plotlabel,color = plt.get_cmap('copper')(float(n)/NT))
            plt.plot(x,Tex, label=plotlabel,color = "green")              
            plt.xlabel(u'$x$', fontsize=26)
            plt.ylabel(u'$T$', fontsize=26, rotation=0)
            plt.title(u'ADRS 1D')
            #plt.legend()
               
        err=np.dot(T-Tex,T-Tex)*dx
        errh1=0
        for j in range (1,NX-1):
            errh1+=dx*(Texx[j]-(T[j+1]-T[j-1])/(2*dx))**2
           
        error[iter]=np.sqrt(err)/NX
        #print('norm error=',error[iter])

        if(abs(time-0.5)<dt*0.5):
            Err_tab1.append(error[iter])

    Err_tab2.append(error[iter])

    
    plt.figure(2)
    plt.plot(np.array(time_tab),rest)

plt.figure(3)
NX_tab=np.array(NX_tab)
Err_tab1=np.array(Err_tab1)
Err_tab2=np.array(Err_tab2)
print(len(NX_tab),len(Err_tab1),len(Err_tab2))

plt.plot(Err_tab1,NX_tab,label="0.5 sec")
plt.plot(Err_tab2,NX_tab,label="1 sec")
plt.ylabel(u'$Nx$', fontsize=14)
plt.xlabel(u'$LÂ² Error$', fontsize=14, rotation=90)
plt.title(u'Error at 2 different times for different meshes')
plt.legend()


# plt.figure(3)
# plt.plot(x,Tex, label=plotlabel,color = plt.get_cmap('copper')(float(n)/NT))
