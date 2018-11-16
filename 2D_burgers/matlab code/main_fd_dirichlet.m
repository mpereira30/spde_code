
clc
clear
close all

T = 2;
dt = 0.01;
a = 2; % a = b case
N = round(T/dt);
fprintf('Number of timesteps: %d \n', N);
J = 30;
h = a/J;
x =(0:h:a)';
y =(0:h:a)';
nu = 0.1; % viscosity
dbcvalue = 0.0; % Dirichlet B.C. value

%  Set initial and desired temperature distributions of rod:
u_0 = zeros(J-1, J-1);
v_0 = zeros(J-1, J-1);
r_x1 = round([0.15*(J-1), 0.75*(J-1)]);
r_x2 = round([0.25*(J-1), 0.85*(J-1)]);
r_y1 = round([0.45*(J-1), 0.45*(J-1)]);
r_y2 = round([0.55*(J-1), 0.55*(J-1)]);

desired_values = [2.0;-2.0];
for i = 1:length(r_x1)
% for i = 2:2
    u_0(r_y1(i):r_y2(i), r_x1(i):r_x2(i)) = desired_values(i);
%     v_0(r_y1(i):r_y2(i), r_x1(i):r_x2(i)) = desired_values(i);
end

num_diff_scheme = 'c'; % b: backward difference and c: central difference

tic
[ut, vt] = pde_fd_dirichlet(u_0(:), v_0(:), h, N, J, dt, nu, dbcvalue, num_diff_scheme); 
% sigma = 0.1;
% [t,ut] = pde_fd_dirichlet_whitenoise(u0(:), T, a, N, J, sigma); % ut = zeros((J-1)*(J-1),length(t));
toc

[X,Y] = meshgrid(x,y);
colormapvalue = cool;

Z_all = sqrt(ut.^2 + vt.^2);
cmin = min(min(Z_all));
cmax = max(max(Z_all));

% plot starting profile
figure(1)
T = 1;
umat = vec2mat(ut(:,T), (J-1))';
umat = padarray(umat, [1 1], dbcvalue, 'both');
vmat = vec2mat(vt(:,T), (J-1))';
vmat = padarray(vmat, [1 1], dbcvalue, 'both');
Z = sqrt(umat.^2 + vmat.^2);
contourf(X,Y,Z);
caxis([cmin cmax]);
colorbar;
hold on;
quiver(X,Y,umat,vmat, 0.5, 'w', 'LineWidth', 0.5);
hold off;
title('start');

% plot midway
figure(2)
T = round((N+1)/2);
umat = vec2mat(ut(:,T), (J-1))';
umat = padarray(umat, [1 1], dbcvalue, 'both');
vmat = vec2mat(vt(:,T), (J-1))';
vmat = padarray(vmat, [1 1], dbcvalue, 'both');
Z = sqrt(umat.^2 + vmat.^2);
contourf(X,Y,Z);
caxis([cmin cmax]);
colorbar;
hold on;
quiver(X,Y,umat,vmat, 0.5, 'w', 'LineWidth', 0.5);
hold off;
title('mid-way');

% plot midway
figure(3)
T = N+1;
umat = vec2mat(ut(:,T), (J-1))';
umat = padarray(umat, [1 1], dbcvalue, 'both');
vmat = vec2mat(vt(:,T), (J-1))';
vmat = padarray(vmat, [1 1], dbcvalue, 'both');
Z = sqrt(umat.^2 + vmat.^2);
contourf(X,Y,Z);
caxis([cmin cmax]);
colorbar;
hold on;
quiver(X,Y,umat,vmat, 0.5, 'w', 'LineWidth', 0.5);
hold off;
title('end');

fig = figure();
for T = 1:N+1
    umat = vec2mat(ut(:,T), (J-1))';
    umat = padarray(umat, [1 1], dbcvalue, 'both');
    vmat = vec2mat(vt(:,T), (J-1))';
    vmat = padarray(vmat, [1 1], dbcvalue, 'both');
% 	Z = sqrt(umat.^2 + vmat.^2);
%     contourf(X,Y,Z);
%     caxis([cmin cmax]);
%     colorbar;
%     hold on;
    quiver(X,Y,umat,vmat, 0.5, 'k', 'LineWidth', 0.5);
%     hold off;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    pause(0.01);
    clf(fig);
end