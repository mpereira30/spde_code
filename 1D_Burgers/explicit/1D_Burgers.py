from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import time, sys

a = 2 # x lies in range (0,a) = (0,2)
nx = 64  # As we increase from 41 to 81, the shape of the initial wave remain conserved more (rectangular instead of gaussian-like)
dx = a / (nx-1)
nt = 500    #nt is the number of timesteps we want to calculate
nu = 0.1 # value of viscosity

periodic = False # else Dirichlet B.C.s

# dt = 0.0095  #dt is the amount of time each timestep covers (delta t)

# Smartly set the time discretization to allow for arbitrary spatial discretization:
sigma = 0.2
# dt = sigma * dx**2 / nu #dt is defined using sigma 
dt = 0.01
print("\ntime discretization is:", dt)
print("end time:", dt*nt)

# Set initial profile:
# u = np.ones(nx)      
u = np.zeros(nx)
u[int(.25 / dx):int(0.75 / dx + 1)] = 2.0
u[int(1.25 / dx):int(1.75 / dx + 1)] = -2.0 

un = np.ones(nx) #initialize a temporary array (store profile at current timestep)

all_u = np.zeros((nt,nx))
for n in range(nt):  # loop for values of n from 0 to nt, so it will run nt times
	un = u.copy() # copy the existing values of u into un
	for i in range(1, nx-1): 
		# if n==0 and i==1:
		# 	print("\nUsing backward difference for advection term!")			
		# u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1]) + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])

		if n==0 and i==1:
			print("\nUsing central difference for advection term!")			
		u[i] = un[i] - un[i] * (0.5*dt / dx) * (un[i+1] - un[i-1]) + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])

	# For periodic B.C.s:
	if(periodic):
		u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 * (un[1] - 2 * un[0] + un[-2])
		u[-1] = u[0]

	all_u[n,:] = u

plt.figure()
ax = plt.gca() 
for i in range(nt):

	sys.stdout.write("step: %d/%d \r" % (i, nt))
	sys.stdout.flush()	

	plt.cla()
	ax.plot(np.linspace(0, 2, nx), all_u[i])
	plt.ylim((-2.5, 2.5))
	plt.title('1D Burgers')
	plt.pause(0.01)

print("Sim complete!")
plt.show()        