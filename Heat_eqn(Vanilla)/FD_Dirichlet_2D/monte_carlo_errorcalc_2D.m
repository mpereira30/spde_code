clc 
clear 

% Based on equation 10.32 in textbook:

sigma = 1;
a = 1;
epsilon = 1;
num_indxs = 1e4;

T = 0.1; % time-instant

u = 0;

for j1 = 1:num_indxs
    for j2 = 1:num_indxs
        lambda_j1_j2 = epsilon*(pi*j1/a)^2 + epsilon*(pi*j2/a)^2;
        u = u + ( (0.5/lambda_j1_j2) * (1 - exp(-2*T*lambda_j1_j2)) );
    end
end

u_theoretical = (sigma^2) * u

M = 1e3; % number of monte carlo samples
J = 4*2.^[2,3,4]'; % spatial-discretization
N = 0.25*J.^2; % time-discretization

for i=1:length(N),
    v(i)=l2_sq_mct(T,a,N(i),J(i),M,sigma);
    disp('index, N, J :')
    disp([i, N(i), J(i)])
    disp([v(i)])
end;


