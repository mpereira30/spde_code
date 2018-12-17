clear 
clc
close all

% set parameters:
T       = 4.0; % total sim time in seconds
dt      = 0.01; 
a       = 2; % rod length 
N       = round(T/dt); % total sim timesteps
nu      = 0.1; % viscosity of medium 
sigma   = 0.05; % space-time noise standard deviation

% spatial discretization: 
% j = 0, 1, 2, ......, J for x = 0, ......., a and x_j = j*h = j*(a/J)
J  = 128; % Total number of indices = J+1 to include j=0 
h  = a/J; % spatial discretization
x  = (0:h:a)';

% set initial profile of velocity:
u0 = zeros(length(x),1);
% u0(round(0.25 / h):round(0.75 / h),1) = 2.0 ; % using the fact that x_j = j*h, therefore, index j = x_j/h
% u0(round(1.25 / h):round(1.75 / h),1) = -2.0; 
range1 = round(0.2*(J+1)):round(0.4*(J+1));
range2 = round(0.6*(J+1)):round(0.8*(J+1));
u0(range1,1) = 0.0;
u0(range2,1) = 0.0;

% Pick numerical method and boundary condition:
method      = 's'; % Different methods are: e (explicit) and s (semi-implicit)
diff_scheme = 'c'; % differentiation scheme for advection: central - c, backward - b
add_noise   = 1;   % 1 - yes, 0 - no

dbc_val_zero = 1.0; 
dbc_val_J    = 1.0;
dbc_val      = [dbc_val_zero, dbc_val_J]; % velocity at boundaries (for Dirichlet B.C.s) ]
u0(1,1)      = dbc_val_zero; % enforce B.C.s at initial time
u0(end,1)    = dbc_val_J; % enforce B.C.s at initial time

tic
ut = pde_fd(u0, dt, h, a, N, J, method, nu, dbc_val, diff_scheme, add_noise, sigma);
toc

figure()
plot(x,ut(:,1));
ylim([-2.0 2.0])
title('start');

figure()
plot(x, ut(:,round(0.25*(N+1))) );
ylim([-2.0 2.0])
title('quarter-way');

figure()
plot(x, ut(:,round(0.5*(N+1))) );
ylim([-2.0 2.0])
title('mid-way');

figure()
plot(x,ut(:,end));
ylim([-2.0 2.0])
title('end');

figure()
t= [0:dt:T];
surf(x,t,ut');

% fig = figure();
% for i = 1:N+1
%     i
%     plot(x, ut(:,i));
%     ylim manual
%     ylim([-2.0 2.0])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.01);
%     clf(fig);
% end

figure()
y = (dt:dt:dt*(N+1));
[X,Y] = meshgrid(x,y);
contourf(X,Y,ut')
colorbar;
caxis([-0.0181 2.6675]);