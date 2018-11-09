from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import time, sys

a = 2 # x lies in range (0,a) = (0,2)
nx = 41  # As we increase from 41 to 81, the shape of the initial wave remain conserved more (rectangular instead of gaussian-like)
dx = a / (nx-1)
nt = 100    #nt is the number of timesteps we want to calculate
dt = .025  #dt is the amount of time each timestep covers (delta t)

go_right = False
if(go_right):
	c = 1.0      # wavespeed
else:
	c = -1.0

# Set initial profile:
u = np.ones(nx)      
u[int(.75 / dx):int(1.25 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s

# plt.figure(1)
# plt.plot(np.linspace(0, 2, nx), u);
# plt.title('initial profile')

un = np.ones(nx) #initialize a temporary array (store profile at current timestep)

all_u = np.zeros((nt,nx))
for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
	un = u.copy() ##copy the existing values of u into un
	
	if(go_right):	# To go right use backward difference in space:

		for i in range(1, nx-1): # this maintains boundary conditions at 1 
			u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])

	else:	# To go left use forward difference in space:

		for i in range(1, nx-1): ## you can try commenting this line and...
			u[i] = un[i] - c * dt / dx * (un[i+1] - un[i])	
	
	all_u[n,:] = u

plt.figure(2)
for i in range(nt):
	plt.cla()
	plt.plot(np.linspace(0, 2, nx), all_u[i])
	plt.ylim((0.95,2.2))
	plt.pause(0.025)

# plt.figure(3)
# plt.plot(np.linspace(0, 2, nx), u);
# plt.title('final profile')

plt.show()        