from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# variable declarations
nx = 64
ny = 64
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx
# dt = 0.01
print("time discretization is:", dt)

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx)) 
v = np.ones((ny, nx))
un = np.ones((ny, nx))
vn = np.ones((ny, nx))

# Assign initial conditions
# set hat function I.C.s : u(.5<=x<=1 && .5<=y<=1 )= v(.5<=x<=1 && .5<=y<=1 ) = 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 
v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

all_us = np.zeros((nt+1,ny,nx))
all_vs = np.zeros((nt+1,ny,nx))
all_us[0,:,:] = u
all_vs[0,:,:] = v
row, col = u.shape

for n in range(nt): 
    un = u.copy()
    vn = v.copy()
    u[1:, 1:] = (un[1:, 1:] - 
                 (un[1:, 1:] * c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                  vn[1:, 1:] * c * dt / dy * (un[1:, 1:] - un[:-1, 1:]))
    v[1:, 1:] = (vn[1:, 1:] -
                 (un[1:, 1:] * c * dt / dx * (vn[1:, 1:] - vn[1:, :-1])) -
                 vn[1:, 1:] * c * dt / dy * (vn[1:, 1:] - vn[:-1, 1:]))
    
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    
    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

    all_us[n+1,:,:] = u
    all_vs[n+1,:,:] = v

X, Y = np.meshgrid(x, y)   
fig = plt.figure(figsize=(11, 7), dpi=100)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.set_xlabel('$x$')
# ax.set_ylabel('$y$')
for n in range(nt):
  plt.cla()
  curr_u = all_us[n,:,:]
  curr_v = all_vs[n,:,:]
  # surf2 = ax.plot_surface(X, Y, curr_u[:], cmap=cm.viridis)
  Q = plt.quiver(X,Y,curr_u, curr_v, units='width')
  plt.pause(0.01)

plt.show()