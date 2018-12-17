clc
clear 
close all

global J;
global N;
global a;
global mu;
global sig;
global T;
global dt;
global sigma;
global h_d;
global rollouts;
global rho;
global scale_factor;
global range;
global epsilon;

% Set parameters:
Tf = 2.5; % It takes around 5 seconds with axon length of 10 for potential to reach the other end
dt = 0.01;
T = round(Tf/dt); % Horizon length
a = 10; % axon length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
iters = 100; % number of PI iterations
rollouts = 100; % number of rollouts for sampling 

rho = 5000;
sigma = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor = 1000;
epsilon = 1; % thermal diffusivity of the rod

range = round(0.7*J):round(0.99*J);

mu = [0.2 0.3 0.4 0.5 0.6 0.7 0.8].*a; % mean location of actuator in space
%mu = [0.2 0.4 0.6 0.8]; % mean location of actuator in space
% mu = [0.7 0.5 0.9].*a; % mean location of actuator in space
sig_val = (0.1*a) * ones(length(mu));
sig = sig_val.^2; % standard deviation of actuation in space

N = length(mu); % number of actuators

h0 = 1./(1+exp(-(2-x)/sqrt(2)));
% h_d(range,1) = 1; % force end of axon to reach potential faster
h_d(range,1) = 0; % suppress axon potential at the end

% Initialize control sequences:
% U = randn(T,N);
U = zeros(T,N);

% First compute matricies curly-M, M and the curly-v-tilde vector offline:
curly_M = compute_curly_M1();
M = compute_capital_M();
curly_v_tilde = compute_curly_v_tilde1();

cost = zeros(iters,1);

f = @(h) h.*(1-h).*(h+0.5);

for iter = 1:iters

    % Generate rollouts and noise samples:
    h_samples = zeros(rollouts,T+1,J+1); % trajectories of spatial points
    xi_samples = zeros(rollouts,T,J+1); % corresponding noise trajectories 
    
    for r = 1:rollouts
        [h_traj, xi_traj] = generate_rollouts1(h0, U, curly_v_tilde, 0, f);
        h_samples(r,:,:) = h_traj;
        xi_samples(r,:,:) = xi_traj;
    end
    
    % Cost computation and control update:
    U_new = PI_control1(h_samples, xi_samples, U, curly_M, M);

    %% Noise-free rollout and noise-free cost computation:
    [h_traj, xi_traj] = generate_rollouts1(h0, U_new, curly_v_tilde, 1, f);

    for t = 1:T
        cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(t,range)' - h_d(range,1))' * (h_traj(t,range)' - h_d(range,1));
    end
    cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(T+1,range)' - h_d(range,1))' * (h_traj(T+1,range)' - h_d(range,1));

    fprintf('Iteration number: %d, Cost: %f \n', iter, cost(iter,1));
    
    U = U_new;
end

figure()
plot(cost)


%%

figure()
plot(x,h_traj(1,:));
title('start');

figure()
plot(x, h_traj(round(0.25*T),:) );
title('quarter-way');

figure()
plot(x, h_traj(round(0.5*T),:) );
title('mid-way');

figure()
plot(x,h_traj(end,:));
hold on
plot(x, h_d, 'g')
title('end');
% 
fig = figure();
for i = 1:T
    plot(x,h_d,'-r');
    hold on;
    plot(x,h_traj(i,:),'-b');
    ylim manual
    ylim([-0.125 1.25])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    pause(0.02);
    clf(fig);
end

figure()
surf(x,(1:1:T+1),h_traj)

figure()
plot(U_new)
