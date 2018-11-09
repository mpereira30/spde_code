from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys
import matplotlib.colors as colors

###variable declarations
nx = 30
ny = 30
nt = 200
c  = 1
a  = 2

# We have nx number of spatial points in the x - direction. However, the indexing goes from 0 to nx -1 
dx = a / (nx - 1) # for domain (0,a), delta x = a/nx-1. Therefore, x_n = n * delta_x = n * a /(nx-1), where n = 0, 1, 2, ...., nx -1
dy = a / (ny - 1)
sigma = .2

nu = 0.1
# dt = sigma * dx * dy / nu
dt = 0.01
print("Parameters used:")
print("nu (viscoity):", nu)
print("time discretization:",dt)
print("r=", nu*dt/(dx**2)) 
print("")

# r is directly proportional to nx or ny

x = np.linspace(0, a, nx)
y = np.linspace(0, a, ny)

dbcvalue = 0.0
u = dbcvalue * np.ones((ny, nx))  # create a 1xn vector of 1's
v = dbcvalue * np.ones((ny, nx))
un = dbcvalue * np.ones((ny, nx)) 
vn = dbcvalue * np.ones((ny, nx))
comb = np.ones((ny, nx))

# Numerical diff for advection term:
use_central_diff  = True
use_backward_diff = False

# Assign initial conditions

# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.75 / dy):int(1.25 / dy + 1),int(.25 / dx):int(0.75 / dx + 1)] = 2 
u[int(.75 / dy):int(1.25 / dy + 1),int(1.25 / dx):int(1.75 / dx + 1)] = -2.0 
# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
# v[int(.25 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

# for i in range(nx):
# 	for j in range(ny):
# 		if (y[j]>=0.5 and y[j]<=1.0) and (x[i]>=0.5 and x[i]<=1.0):
# 			u[i,j] = 0.0
# 			v[i,j] = 1.0
# 		else:
# 			u[i,j] = 1.0
# 			v[i,j] = 0.0

all_us = np.zeros((nt+1,ny,nx))
all_vs = np.zeros((nt+1,ny,nx))
all_us[0,:,:] = u
all_vs[0,:,:] = v

for n in range(nt): ##loop across number of time steps
	un = u.copy()
	vn = v.copy()

	if use_backward_diff:
		if n == 0:
			print("Using backward difference for advection")
		u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
		                 dt / dx * un[1:-1, 1:-1] * 
		                 (un[1:-1, 1:-1] - un[1:-1, 0:-2]) - 
		                 dt / dy * vn[1:-1, 1:-1] * 
		                 (un[1:-1, 1:-1] - un[0:-2, 1:-1]) + 
		                 nu * dt / dx**2 * 
		                 (un[1:-1,2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + 
		                 nu * dt / dy**2 * 
		                 (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

		v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - 
		                 dt / dx * un[1:-1, 1:-1] *
		                 (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
		                 dt / dy * vn[1:-1, 1:-1] * 
		                (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) + 
		                 nu * dt / dx**2 * 
		                 (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
		                 nu * dt / dy**2 *
		                 (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))
	elif use_central_diff:
		if n == 0:
			print("Using backward difference for advection")		
		u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
		                 (0.5*dt / dx) * un[1:-1, 1:-1] * 
		                 (un[1:-1, 2:] - un[1:-1, 0:-2]) - 
		                 (0.5*dt / dy) * vn[1:-1, 1:-1] * 
		                 (un[2:, 1:-1] - un[0:-2, 1:-1]) + 
		                 nu * dt / dx**2 * 
		                 (un[1:-1,2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + 
		                 nu * dt / dy**2 * 
		                 (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

		v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - 
		                 (0.5*dt / dx) * un[1:-1, 1:-1] *
		                 (vn[1:-1, 2:] - vn[1:-1, 0:-2]) -
		                 (0.5*dt / dy) * vn[1:-1, 1:-1] * 
		                (vn[2:, 1:-1] - vn[0:-2, 1:-1]) + 
		                 nu * dt / dx**2 * 
		                 (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
		                 nu * dt / dy**2 *
		                 (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))		
	 
	u[0, :] = dbcvalue
	u[-1, :] = dbcvalue
	u[:, 0] = dbcvalue
	u[:, -1] = dbcvalue

	v[0, :] = dbcvalue
	v[-1, :] = dbcvalue
	v[:, 0] = dbcvalue
	v[:, -1] = dbcvalue

	all_us[n+1,:,:] = u
	all_vs[n+1,:,:] = v

fig = plt.figure(figsize=(11, 7), dpi=100)
# ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
# ax.set_xlabel('$x$')
# ax.set_ylabel('$y$');

all_Z = np.sqrt(all_us**2 + all_vs**2)
vmax = np.amax(all_Z)
vmin = np.amin(all_Z)

# print(vmax, vmin)

# print(np.amax(all_us), np.amax(all_vs), np.amin(all_us), np.amin(all_vs))

# fig = plt.figure(figsize=(11, 7), dpi=100)
# ax = fig.gca(projection='3d')

for n in range(nt):
	# print(n)
	plt.cla()
	curr_u = all_us[n,:,:]
	curr_v = all_vs[n,:,:]

	# ax.clear()
	# ax.plot_surface(X, Y, curr_u[:], cmap=cm.viridis, rstride=1, cstride=1)
	# ax.plot_surface(X, Y, curr_v[:], cmap=cm.viridis, rstride=1, cstride=1)	
	# ax.set_zlim3d(0.95,2.05)

	sys.stdout.write("t: %f/%f \r" %(n*dt,nt*dt))
	sys.stdout.flush()
	Z = np.sqrt(curr_u**2 + curr_v**2)
	# plt.contourf(X, Y, Z, cmap='Blues', vmin=vmin, vmax=vmax)
	# if n==0:
	# 	plt.colorbar();
	# Q = plt.quiver(X,Y,curr_u, curr_v, units='xy', angles='xy', scale_units='xy', scale=10, color='black')		
	Q = plt.quiver(X,Y,curr_u, curr_v, Z, cmap='autumn', clim= [vmin, vmax], units='xy', angles='xy', scale_units='xy', scale=15)
	if n==0:
		plt.colorbar(Q,extend='max');		

	plt.pause(0.001)


plt.show()