
clc
clear
close all

T = 5;
dt = 0.01;
a = [5 5]; % a = b case
N = round(T/dt)
J = [50 50];

x=(0 : a(1)/J(1) : a(1))';  % hx = ax/J, x contains (J+1) elements where x = a(1) is (J+1)th element
y=(0 : a(2)/J(2) : a(2))';  % hy = ay/J, y contains (J+1) elements where y = a(2) is (J+1)th element

init_temp = 5;
u0 = init_temp*ones((J(1)+1)*(J(2)+1),1);
Len = length(u0);
u0(1 : round(0.4*Len))=0;
u0(round(0.6*Len) : end)=0;
u0 = vec2mat(u0,(J(1)+1)); 
u0 = u0';

tic
alpha = 0.1; % called regularity parameter. Used to compute the 2D eigenvalues exp(-alpha * lambda_j1_j2)
[t,u,ut] = spde_twod_Gal(u0,T,a,N,J,alpha); 
toc

figure(1)
surf(x,y,ut(:,:,1));
title('at t=0');
xlabel('x');
ylabel('y');

figure(2)
surf(x,y,ut(:,:,end));
title('at t=T');
xlabel('x');
ylabel('y');

% fig = figure(3);
% for i = 1:length(t)
%     surf(x,y,ut(:,:,i));
%     zlim manual
%     zlim([0 1e-3])    
% %     zlim([0 init_temp])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.01);
%     clf(fig);
% end