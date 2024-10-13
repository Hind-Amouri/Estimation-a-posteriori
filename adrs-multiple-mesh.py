import math
import numpy as np
import matplotlib.pyplot as plt

#u,t = -V u,x + k u,xx  -lamda u + f

# PHYSICAL PARAMETERS
K = 0.1     #Diffusion coefficient
L = 1.0     #Domain size
Time = 20.  #Integration time


V=1
lamda=1

# NUMERICAL PARAMETERS
NX = 10    #Number of grid points
NT = 10000   #Number of time steps max
ifre=1000000  #plot every ifre time iterations
eps=0.001     #relative convergence ratio
niter_refinement=10      #niter different calculations

errorL2=np.zeros((niter_refinement))
errorH1=np.zeros((niter_refinement))
semiH2=np.zeros((niter_refinement))
itertab=np.zeros((niter_refinement))

for iter in range (niter_refinement):
	NX=NX+5

	dx = L/(NX-1)                 #Grid step (space)
	dt = dx**2/(V*dx+4*K+dx**2)   #Grid step (time)  condition CFL de stabilite 10.4.5
	print(dx,dt)
	itertab[iter]=dx

	### MAIN PROGRAM ###

	# Initialisation
	x = np.linspace(0.0,1.0,NX)
	T = np.zeros((NX)) #np.sin(2*np.pi*x)
	F = np.zeros((NX))
	rest = []
	RHS = np.zeros((NX))

	Tex = np.zeros((NX)) #np.sin(2*np.pi*x)
	for j in range (1,NX-1):
		Tex[j] = np.exp(-20*(j*dx-0.5)**2)
	for j in range (1,NX-1):
		Tx=(Tex[j+1]-Tex[j-1])/(2*dx)
		Txx=(Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)
		F[j]=V*Tx-K*Txx+lamda*Tex[j]
		

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


		if (n == 1 ):
			res0=res

		rest.append(res)
	#Plot every ifre time steps
		if (n%ifre == 0 or (res/res0)<eps):
			print(n,res)
			plotlabel = "t = %1.2f" %(n * dt)
			plt.plot(x,T, label=plotlabel,color = plt.get_cmap('copper')(float(n)/NT))
		  

	print(n,res)
	plt.plot(x,T)

	plt.xlabel(u'$x$', fontsize=26)
	plt.ylabel(u'$T$', fontsize=26, rotation=0)
	plt.title(u'ADRS 1D')
	plt.legend()

	plt.figure(2)
	plt.plot(np.log10(rest/rest[0]))


#	errL2=np.sqrt(np.dot(T-Tex,T-Tex))
	errH1h=0
	errL2h=0
	semih2=0
	for j in range (1, NX-1):
		Texx=(Tex[j+1]-Tex[j-1])/(2*dx)
		Tx=(T[j+1]-T[j-1])/(2*dx)
		errL2h+=dx*(T[j]-Tex[j])**2
		errH1h+=dx*(Tx-Texx)**2
		Txx=(Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)
		semih2+=dx*Txx**2

	errorL2[iter]=errL2h
	errorH1[iter]=errL2h+errH1h
	semiH2[iter]=semih2
	
	
	print('norm error L2, H1=',errL2h,errH1h)


plt.figure(3)
#plt.plot(itertab,np.log10(errorL2),itertab,np.log10(errorH1))

plt.plot(itertab,np.sqrt(errorL2)/semiH2,itertab,np.sqrt(errorH1)/semiH2)



plt.show()
