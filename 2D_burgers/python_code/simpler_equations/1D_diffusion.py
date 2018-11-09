from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import time, sys

a = 2 # x lies in range (0,a) = (0,2)
nx = 41  # As we increase from 41 to 81, the shape of the initial wave remain conserved more (rectangular instead of gaussian-like)
dx = a / (nx-1)
nu = 0.3 # value of viscosity
nt = 100    #nt is the number of timesteps we want to calculate
sigma = 0.2 
dt = sigma * dx**2 / nu #dt is defined using sigma 
print("\ntime discretization is:", dt)

# Set initial profile:
u = np.ones(nx)      
u[int(.75 / dx):int(1.25 / dx + 1)] = 2 

un = np.ones(nx) #initialize a temporary array (store profile at current timestep)

all_u = np.zeros((nt,nx))
for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
	un = u.copy() ##copy the existing values of u into un

	for i in range(1, nx-1): 
		u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
	
	all_u[n,:] = u

plt.figure()
for i in range(nt):
	plt.cla()
	plt.plot(np.linspace(0, 2, nx), all_u[i])
	plt.ylim((0.95,2.2))
	plt.pause(0.01)

plt.show()        