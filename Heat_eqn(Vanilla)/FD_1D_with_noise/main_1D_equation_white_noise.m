clear 
clc
close all

% set parameters:
T = 5;
dt = 0.01;
a = 5; % rod length 
N = round(T/dt)

% spatial discretization: 
% j = 0, 1, 2, ......, J for x = 0, ......., a and x_j = j*h = j*(a/J)
J = 256; 
h = a/J;
x = (0:h:a)';

% set initial temperature profile:
% u0=1./(1+exp(-(2-x)/sqrt(2)));
u0 = zeros(length(x),1);
init_temp = 50;
u0(2:round(0.25*J)) = init_temp;

bctype = 'D' % select boundary condition type
tic
sigma = 1.5;
epsilon = 1.0;
% [t, ut] = pde_fd_white_noise(u0,T,a,N,J,bctype,sigma, epsilon);
[t, ut] = pde_fd_white_noise_series(u0,T,a,N,J,bctype,sigma, epsilon); % ONLY DIRICHELT BCs ! 
toc

figure()
plot(x,ut(:,1));
title('start');

figure()
plot(x, ut(:,round(0.25*length(t))) );
title('quarter-way');

figure()
plot(x, ut(:,round(0.5*length(t))) );
title('mid-way');

figure()
plot(x,ut(:,end));
title('end');


% fig = figure();
% for i = 1:length(t)
%     plot(x, ut(:,i));
%     ylim manual
% %     zlim([0 1e-3])    
%     ylim([0 init_temp])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.01);
%     clf(fig);
% end