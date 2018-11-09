
clc
clear
close all

T = 5;
dt = 0.01;
a = 5; % a = b case
N = round(T/dt)
J = 128;
h = a/J;
x =(0:h:a)';
y =(0:h:a)';

% heat a rectangular portion of the plate:
init_temp = 1;
u0 = zeros((J-1),(J-1));
u0( 1:round(0.2*(J-1)), 1:round(0.2*(J-1)) ) = init_temp;
u0( round(0.8*(J-1)):end, round(0.8*(J-1)):end ) = init_temp;
u0( round(0.8*(J-1)):end, 1:round(0.2*(J-1)) ) = init_temp;
u0( 1:round(0.2*(J-1)), round(0.8*(J-1)):end ) = init_temp;

% u0( 1:round(0.2*(J-1)), round(0.4*(J-1)):round(0.6*(J-1)) ) = init_temp;
% u0( round(0.8*(J-1)):end, round(0.4*(J-1)):round(0.6*(J-1)) ) = - init_temp;

% Gaussian temp distribution
% h = fspecial('gaussian', [(J-1) (J-1)], 10);
% u0 = h(:);

tic
% [t,ut] = pde_fd_dirichlet(u0(:), T, a, N, J); % ut = zeros((J-1)*(J-1),length(t)); 
sigma = 0.1;
[t,ut] = pde_fd_dirichlet_whitenoise(u0(:), T, a, N, J, sigma); % ut = zeros((J-1)*(J-1),length(t));
toc

figure()
Z = vec2mat(ut(:,1), (J-1))';
Z = padarray(Z, [1 1], 0, 'both');
surf(x,y,Z);
title('start');

figure()
Z = vec2mat(ut(:,round(0.5*length(t))), (J-1))';
Z = padarray(Z, [1 1], 0, 'both');
surf(x,y,Z);
title('half-way');

figure()
Z = vec2mat(ut(:,end), (J-1))';
Z = padarray(Z, [1 1], 0, 'both');
surf(x,y,Z);
title('end');
% % 
% fig = figure();
% for i = 1:length(t)
%     Z = vec2mat(ut(:,i), (J-1))';
%     Z = padarray(Z, [1 1], 0, 'both');
%     surf(x,y,Z);
%     zlim manual
% %     zlim([0 1e-3])    
%     zlim([0 init_temp])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.01);
%     clf(fig);
% end