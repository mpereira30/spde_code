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
global terminal_only;

terminal_only = 0; % Set 1 for considering only terminal state cost

% Set parameters:
dt = 0.01;

% Total simulation time and timesteps:
Tsim = 1.5; % in seconds
Tsim_steps = round(Tsim/dt);

% MPC final time horizon and timesteps
Tf = 0.1; % in seconds
T = round(Tf/dt); 

a = 10; % axon length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
iters = 5; % number of PI iterations
rollouts = 100; % number of rollouts for sampling 

rho = 10000;
sigma = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor = 1000;
epsilon = 1; % thermal diffusivity of the rod

range = round(0.7*J):round(0.99*J);

mu = [0.2 0.3 0.4 0.5 0.6 0.7 0.8].*a; % mean location of actuator in space
sig_val = (0.1*a) * ones(length(mu));
sig = sig_val.^2; % standard deviation of actuation in space

N = length(mu); % number of actuators

h0 = 1./(1+exp(-(2-x)/sqrt(2)));

h_d = zeros(size(x));
h_d(range,1) = 1; % force end of axon to reach potential faster
% h_d(range,1) = 0; % suppress axon potential at the end

% Initialize control sequences:
U = randn(T,N);
% U = zeros(T,N);

% First compute matricies curly-M, M and the curly-v-tilde vector offline:
curly_M = compute_curly_M1();
M = compute_capital_M();
curly_v_tilde = compute_curly_v_tilde1();

cost = zeros(Tsim_steps,1);

f = @(h) h.*(1-h).*(h+0.5);

% Container to store the trajectory of the actual system:
h_actual = zeros(Tsim_steps+1, J+1);
h_actual(1,:) = h0; % store the initial temperature profile

for cur_step = 1:Tsim_steps

    fprintf('current timestep: %d of %d, ',cur_step, Tsim_steps );
    
    %% MPC iterations:
    for iter = 1:iters  
        
        % Generate rollouts and noise samples:
        h_samples = zeros(rollouts,T+1,J+1); % trajectories of spatial points
        xi_samples = zeros(rollouts,T,J-1); % corresponding noise trajectories 

        for r = 1:rollouts
            [h_traj, xi_traj] = generate_rollouts1(h0, U, curly_v_tilde, 0, f);
            h_samples(r,:,:) = h_traj;
            xi_samples(r,:,:) = xi_traj;
        end

        % Cost computation and control update:
        U_new = PI_control1(h_samples, xi_samples, U, curly_M, M);
        U = U_new;
       
    end
    
    % Play on NOISY DYNAMICS:
    [h_traj, xi_traj] = generate_rollouts1(h0, U_new, curly_v_tilde, 1, f);
    
    % Take the propagated state after application of 1st control of U_new:
    h0 = ( h_traj(2,:) )';
    h_actual(cur_step+1, :) = h0';
    
    % Set the initial control sequence for the next time step:
    U = [U_new(2:end,:); U_new(end,:)];
    
    cost(cur_step) = scale_factor * (h0(range,1) - h_d(range,1))' * (h0(range,1) - h_d(range,1));
    fprintf('cost: %f\n', cost(cur_step));
end

figure()
plot(cost)


%%

figure()
plot(x,h_actual(1,:));
title('start');

figure()
plot(x, h_actual(round(0.25*Tsim_steps),:) );
title('quarter-way');

figure()
plot(x, h_actual(round(0.5*Tsim_steps),:) );
title('mid-way');

figure()
plot(x,h_actual(end,:));
hold on
plot(x, h_d, 'g')
title('end');
% 
% fig = figure();
% for i = 1:Tsim_steps+1
%     plot(x,h_d,'-r');
%     hold on;
%     plot(x,h_actual(i,:),'-b');
%     ylim manual
%     ylim([-0.125 1.25])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.02);
%     clf(fig);
% end

figure()
surf((1:1:Tsim_steps+1),x,h_actual')

figure()
plot(U_new)

figure()
contourf(h_actual)