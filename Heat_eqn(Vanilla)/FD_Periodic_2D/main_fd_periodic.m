
clc
clear
close all

T = 5;
dt = 0.01;
a = 5; % a = b case
N = round(T/dt)
J = 100;
h = a/J;
x =(0:h:a)';
y =(0:h:a)';

% heat a rectangular portion of the plate:
init_temp = 5;
u0 = zeros(J,J);
% u0( 1:round(0.2*J), 1:round(0.2*J) ) = init_temp;
% u0( round(0.8*J):end, round(0.8*J):end ) = init_temp;
% u0( round(0.8*J):end, 1:round(0.2*J) ) = init_temp;
% u0( 1:round(0.2*J), round(0.8*J):end ) = init_temp;

u0( 1:round(0.2*J), round(0.4*J):round(0.6*J) ) = init_temp;
u0( round(0.8*J):end, round(0.4*J):round(0.6*J) ) = - init_temp;

% u0( 1:round(0.2*J), round(0.4*J):round(0.6*J) ) = init_temp;
% u0( round(0.8*J):end, round(0.4*J):round(0.6*J) ) = init_temp;

% Gaussian temp distribution
% h = fspecial('gaussian', [(J+1) (J+1)], 10);
% h = h(1:J,1:J);
% u0 = h(:);

tic
[t,ut] = pde_fd_periodic(u0(:), T, a, N, J); % ut = zeros(J*J,length(t)); 
toc

figure()
Z = vec2mat(ut(:,1), J)';
Z = [Z, Z(:,1)]; 
Z = [Z; Z(1,:)];
surf(x,y,Z);
title('start');

figure()
Z = vec2mat(ut(:,round(0.5*length(t))), J)';
Z = [Z, Z(:,1)]; 
Z = [Z; Z(1,:)];
surf(x,y,Z);
title('half-way');

figure()
Z = vec2mat(ut(:,end), J)';
Z = [Z, Z(:,1)]; 
Z = [Z; Z(1,:)];
surf(x,y,Z);
title('end');

fig = figure(3);
for i = 1:length(t)
    Z = vec2mat(ut(:,i), J)';
    Z = [Z, Z(:,1)]; 
    Z = [Z; Z(1,:)];
    surf(x,y,Z);
    zlim manual
%     zlim([0 1e-3])    
    zlim([-init_temp init_temp])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    pause(0.01);
    clf(fig);
end