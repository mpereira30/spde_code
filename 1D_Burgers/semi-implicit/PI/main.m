clc
clear 
close all

% TODOs :
% (1.) What happens when sigma changes? Should we instead multiply the Xi directly by sigma instead of multiplying the dW by sigma?
% (2.) In computation for zeta_1, sqrt(dt) or dt ??
% (3.) In control update equation, sqrt(dt) or dt ??

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
Tf = 1.0; % final time (in seconds)
dt = 0.01;
T = round(Tf/dt); % Horizon length
a = 2; % rod length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
iters = 1000; % number of PI iterations
rollouts = 200; % number of rollouts for sampling 

rho = 1; % temperature parameter
sigma = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor = 100;
epsilon = 1; % thermal diffusivity of the rod

%  Set initial and desired temperature distributions of rod:
init_temp = 0;
desired_temp = 5;

h_d = zeros(J+1,1);

range1 = round(0.18*J):round(0.22*J);
range2 = round(0.48*J):round(0.52*J);
range3 = round(0.78*J):round(0.82*J);
range = [range1, range2, range3];
h_d(range1,1) = +desired_temp;
h_d(range2,1) = +desired_temp * 0.5;
h_d(range3,1) = +desired_temp;

mu = [0.2 0.5 0.8].*a;
sig_val = (0.1*a) * ones(length(mu));
sig = sig_val.^2; % standard deviation of actuation in space

N = length(mu); % number of actuators

% Set the initial temperature profile of the rod:
h0 = zeros(J+1,1);

% Initialize control sequences:
U = randn(T,N);

% First compute matricies curly-M, M and the curly-v-tilde vector offline:
curly_M = compute_curly_M();
M = compute_capital_M();
curly_v_tilde = compute_curly_v_tilde();

cost = zeros(iters,1);

for iter = 1:iters

    % Generate rollouts and noise samples:
    h_samples = zeros(rollouts,T+1,J+1); % trajectories of spatial points
    xi_samples = zeros(rollouts,T,J-1); % corresponding noise trajectories 
    
    for r = 1:rollouts
        [h_traj, xi_traj] = generate_rollouts(h0, U, curly_v_tilde, 0);
        h_samples(r,:,:) = h_traj;
        xi_samples(r,:,:) = xi_traj;
    end
    
    % Cost computation and control update:
    U_new = PI_control(h_samples, xi_samples, U, curly_M, M);

    %% Noise-free rollout and noise-free cost computation:
%     [h_traj, xi_traj] = generate_rollouts(h0, U_new, curly_v_tilde, 1);
    
    % With disturbances in the actual dynamics:
    [h_traj, xi_traj] = generate_rollouts(h0, U_new, curly_v_tilde, 0);    
    
    if (terminal_only == 0) % consider running cost as well
        for t = 1:T
            cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(t,range)' - h_d(range,1))' * (h_traj(t,range)' - h_d(range,1));
        end
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
title('end');
% 
% fig = figure();
% for i = 1:T
%     plot(x,h_d, '-r');
%     hold on;
%     plot(x,h_traj(i,:), '-b');
%     ylim manual
% %     zlim([0 1e-3])    
% %     ylim([-desired_temp*1.25 desired_temp*1.25])
%     ylim([0 desired_temp*1.25])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.03);
%     clf(fig);
% end

figure()
surf(x,(1:1:T+1),h_traj)

 
figure()
plot(U_new)

mean_h = mean(h_traj, 1);
var_h = var(h_traj,1); 
figure()
plot(x,h_d, '-r');
hold on; 
plot(x, mean_h,'-b');
hold on; 
plot(x, mean_h + 2.*sqrt(var_h),'--k');
hold on; 
plot(x, mean_h - 2.*sqrt(var_h),'--k');