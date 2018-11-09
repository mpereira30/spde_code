from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
# dt = sigma * dx
dt = 0.01
print("time discretization is:", dt)

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx)) 
un = np.ones((ny, nx))

# Assign initial conditions
# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 

# Plot Initial Condition:
# fig = plt.figure(figsize=(11, 7), dpi=100)
# ax = fig.gca(projection='3d') # you need Axes 3D for this !                      
# X, Y = np.meshgrid(x, y)                            
# surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
# plot_surface  equivalent to the regular plot command, but it takes a grid of X and Y values for the data point positions

all_us = np.zeros((nt+1,ny,nx))
all_us[0,:,:] = u
row, col = u.shape

# Naive way (nested for-loops):
# for n in range(nt): 
#     un = u.copy()
#     for j in range(1, row):
#         for i in range(1, col):
#             u[j, i] = (un[j, i] - (c * dt / dx * (un[j, i] - un[j, i - 1])) -
#                                   (c * dt / dy * (un[j, i] - un[j - 1, i])))
#             u[0, :] = 1
#             u[-1, :] = 1
#             u[:, 0] = 1
#             u[:, -1] = 1

#     all_us[n+1,:,:] = u

# Efficient way:
for n in range(nt): ##loop across number of time steps
    un = u.copy()
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                              (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    all_us[n+1,:,:] = u

# fig = plt.figure(figsize=(11, 7), dpi=100)
X, Y = np.meshgrid(x, y)   
fig = plt.figure()
ax = fig.gca(projection='3d')
plt.xlabel('x')
plt.ylabel('y')
for n in range(nt):
	plt.cla()
	curr_u = all_us[n,:,:]
	surf2 = ax.plot_surface(X, Y, curr_u[:], cmap=cm.viridis)
	plt.pause(0.01)

plt.show()