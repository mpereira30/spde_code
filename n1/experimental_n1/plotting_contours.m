clc
clear
close all

dt = 0.01;
a = 5.0; % axon length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
Tsim_suppress       = 5.0; % in seconds
Tsim_accel          = 1.5; % in seconds
Tsim_steps_accel    = round(Tsim_accel/dt);
Tsim_steps_suppress = round(Tsim_suppress/dt);

load('suppress_data.mat');
load('accel_data');
load('vanilla_nagumo.mat');

h_all = [h_accel; h_suppress; h_traj];
cmin = min(min(h_all));
cmax = max(max(h_all));

figure()
y = (dt:dt:dt*(Tsim_steps_accel+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,h_accel)
colorbar;
caxis([cmin cmax]);
ylabel('time in seconds');
xlabel('spatial position along axon');
title('accel');

figure()
y = (dt:dt:dt*(Tsim_steps_suppress+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,h_suppress)
colorbar;
caxis([cmin cmax]);
title('suppress');
ylabel('time in seconds');
xlabel('spatial position along axon');

figure()
y = (dt:dt:dt*(Tsim_steps_suppress+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,h_traj)
colorbar;
caxis([cmin cmax]);
title('vanilla');
ylabel('time in seconds');
xlabel('spatial position along axon');

% figure()
% subplot(1,3,1)
% y = (dt:dt:dt*(Tsim_steps_suppress+1));
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,h_traj)
% title('vanilla');
% 
% subplot(1,3,2)
% y = (dt:dt:dt*(Tsim_steps_accel+1));
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,h_accel)
% title('accel');
% 
% subplot(1,3,3)
% y = (dt:dt:dt*(Tsim_steps_suppress+1));
% [X,Y] = meshgrid(x,y);
% contourf(X,Y,h_suppress)
% colorbar;
% caxis([cmin cmax]);
% title('suppress');

