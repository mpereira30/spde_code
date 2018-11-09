clc 
clear 

% Based on equation 10.43 in textbook:

sigma = 1;
a = 1;
epsilon = 1;
num_indxs = 1e6;

T = 0.1; % time-instant

u = 0;

for j = 1:num_indxs
    lambda_j = epsilon*(pi*j/a)^2;
    u = u + ( (0.5/lambda_j) * (1 - exp(-2*T*lambda_j)) );
end

u_theoretical = (sigma^2) * u

M = 1e4; % number of monte carlo samples
J = 4*2.^[2,3,4]'; % spatial-discretization
N = 0.25*J.^2; % time-discretization
bctype = 'D';

% disp('white noise with kernel')
% 
% for i=1:length(N),
%     v(i) = l2_sq_mct(T,a,N(i),J(i),M,epsilon,sigma,bctype);
%     disp('index, N, J :')
%     disp([i, N(i), J(i)])
%     disp([v(i)])
% end;

disp('white noise with series')

for i=1:length(N),
    v(i) = l2_sq_mct_series(T,a,N(i),J(i),M,epsilon,sigma,bctype);
    disp('index, N, J :')
    disp([i, N(i), J(i)])
    disp([v(i)])
end;



