clc
clear
close all

dt = 0.01;
a = 10; % axon length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
Tsim = 1.5; % in seconds
Tsim_steps = round(Tsim/dt);

load('suppress_data.mat');
load('accel_data');
load('vanilla_nagumo.mat');

h_all = [h_accel; h_suppress; h_traj];
cmin = min(min(h_all));
cmax = max(max(h_all));

figure()
y = (dt:dt:dt*(Tsim_steps+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,h_accel)
colorbar;
caxis([cmin cmax]);
title('accel');

figure()
y = (dt:dt:dt*(Tsim_steps+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,h_suppress)
colorbar;
caxis([cmin cmax]);
title('suppress');

figure()
y = (dt:dt:dt*(1000+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,h_traj)
colorbar;
caxis([cmin cmax]);
title('vanilla');
