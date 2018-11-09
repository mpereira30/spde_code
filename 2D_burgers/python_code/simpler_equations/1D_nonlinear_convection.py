from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import time, sys

a = 2 # x lies in range (0,a) = (0,2)
nx = 30  # As we increase from 41 to 81, the shape of the initial wave remain conserved more (rectangular instead of gaussian-like)
dx = a / (nx-1)
nt = 100    #nt is the number of timesteps we want to calculate
# dt = .025  #dt is the amount of time each timestep covers (delta t)

dt = 0.001

# Set initial profile:
go_right = True
if(go_right):
	u = np.ones(nx)      
	u[int(.75 / dx):int(1.25 / dx + 1)] = 2 
else:
	u = -1.0 * np.ones(nx)      
	u[int(.75 / dx):int(1.25 / dx + 1)] = -2.0  

un = np.ones(nx) #initialize a temporary array (store profile at current timestep)

all_u = np.zeros((nt,nx))
for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
	un = u.copy() ##copy the existing values of u into un
	if(go_right):
		for i in range(1, nx-1): ## you can try commenting this line and...
			# if n==0 and i==1:
			# 	print("\nUsing backward difference !")		
			# u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1])

			# Trying central difference instead of backward difference:
			# if n==0 and i==1:
			# 	print("\nUsing central difference !")
			# u[i] = un[i] - un[i] * (0.5 * dt / dx) * (un[i+1] - un[i-1])

			# Trying central difference instead of backward difference:
			if n==0 and i==1:
				print("\nUsing forward difference !")
			u[i] = un[i] - un[i] * (dt / dx) * (un[i+1] - un[i])			

	else:
		for i in range(1, nx-1): ## you can try commenting this line and...
			u[i] = un[i] - un[i] * dt / dx * (un[i+1] - un[i])
	
	all_u[n,:] = u

plt.figure()
for i in range(nt):
	plt.cla()
	plt.plot(np.linspace(0, 2, nx), all_u[i])
	plt.ylim((0.95,2.2))
	plt.pause(0.02)

print("sim complete!")
plt.show()        