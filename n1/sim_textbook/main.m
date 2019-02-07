clc
clear 
close all

global J;
global a;
global T;
global dt;
global sigma;
global epsilon;

dt = 0.01;
a = 5.0; 
J = 128;
x = (0:a/J:a)'; 
ell = 1; 
Tf = 5.0;
T = round(Tf/dt);
epsilon = 1; 
rho = 10000;
sigma = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
noise_free = 0;

u0 = 1./(1+exp(-(2-x)/sqrt(2)));
[h_traj, xi_traj] = generate_rollouts( u0, noise_free, @(u) u .* (1-u) .* (u+0.5) );

% figure()
% plot(x,h_traj(1,:));
% title('start');
% 
% figure()
% plot(x, h_traj(round(0.25*T),:) );
% title('quarter-way');
% 
% figure()
% plot(x, h_traj(round(0.5*T),:) );
% title('mid-way');
% 
% figure()
% plot(x,h_traj(end,:));
% title('end');
% 
% fig = figure();
% for i = 1:T
%     plot(x,h_traj(i,:), '-b');
%     ylim manual
%     ylim([0 1])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.0051);
%     clf(fig);
% end

% figure()
% surf(x,(1:1:T+1),h_traj)

figure()
y = (dt:dt:dt*(T+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y, h_traj);
xlabel('x-position')
ylabel('time')