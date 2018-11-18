clc
clear 
close all

global J N a mu sig T dt sigma h_d rollouts rho scale_factor range nu terminal_only 

%-------------- Setup simulation parameters -------------------------------

Tf              = 1.0; % final time (in seconds)
dt              = 0.01;
T               = round(Tf/dt); % Horizon length or total number of timesteps
a               = 2; % rod length 
J               = 128; % number of spatial points
z               = a/J; % spatial discretization
x               = (0:z:a)';
iters           = 100; % number of PI iterations (Open Loop)
rollouts        = 100; % number of rollouts for sampling 
rho             = 500; % Path Integral temperature parameter, High rho => like max fn, Low rho => like averaging the samples
sigma           = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor    = 100; % For scaling certain terms in cost function to increase relative importance
nu              = 0.1; % viscosity of medium
terminal_only   = 1; % Set 1 for considering only terminal state cost
load_U          = 0; % Set to 1 if you want to load previous optimized control sequence

%------- Set initial and desired velocity profiles ------------------------

desired_vel     = 2.0;
h0              = zeros(J+1,1);
h_d             = NaN(J+1,1);

range1          = round(0.18*(J+1)):round(0.22*(J+1));
range2          = round(0.48*(J+1)):round(0.52*(J+1));
range3          = round(0.78*(J+1)):round(0.82*(J+1));

h_d(range1,1)   = desired_vel;
h_d(range2,1)   = desired_vel * 0.5;
h_d(range3,1)   = +desired_vel;

range           = [range1, range2, range3];
% range           = [range1, range2];

%----------- Set the Dirichlet B.C. values at each end --------------------

dbc_val_zero = 1.0; % B.C. at start of spatial domain
dbc_val_J    = 1.0; % B.C. at end of spatial domain
dbc_val      = [dbc_val_zero, dbc_val_J]; 
h0(1,1)      = dbc_val_zero; % enforce B.C.s at initial time
h0(end,1)    = dbc_val_J; % enforce B.C.s at initial time

%--------------------- Setup actuators ------------------------------------

mu           = [0.2 0.5 0.8] .* a; 
sig_val      = (0.05*a) * ones(length(mu));
sig          = sig_val.^2; % standard deviation of actuation in space
N            = length(mu); % number of actuators

%---------- Initialize controls and compute offline matrices -------------- 

if load_U == 1
    fprintf("Loading previously optimized control sequence\n");
    load('optimal_controls.mat')
    U           = U_new;
else
    U           = randn(T,N);
end
curly_M         = compute_curly_M();
M               = compute_capital_M();
curly_v_tilde   = compute_curly_v_tilde();

cost = zeros(iters,1);

% Local copies of global variables to be used in parfor loop:
J_     = J;
a_     = a;
T_     = T;
dt_    = dt;
nu_    = nu;
sigma_ = sigma;

for iter = 1:iters

    % Generate rollouts and noise samples:
    h_samples = zeros(rollouts,T+1,J+1); % trajectories of spatial points
    xi_samples = zeros(rollouts,T,J-1); % corresponding noise trajectories 
    
    noise_free = 0; % add noise for exploration
    parfor r = 1:rollouts
        [h_traj, xi_traj] = generate_rollouts(h0, U, curly_v_tilde, noise_free, J_, nu_, a_, T_, dbc_val, dt_, sigma_);
        h_samples(r,:,:) = h_traj;
        xi_samples(r,:,:) = xi_traj;
    end
    
    % Cost computation and control update:
    U_new = PI_control(h_samples, xi_samples, U, curly_M, M);

    %% Noise-free rollout and noise-free cost computation:
    noise_free = 1;
    [h_traj, xi_traj] = generate_rollouts(h0, U_new, curly_v_tilde, noise_free, J_, nu_, a_, T_, dbc_val, dt_, sigma_);
    
    % With disturbances in the actual dynamics:
%     [h_traj, xi_traj] = generate_rollouts(h0, U_new, curly_v_tilde, 0);    
    
    if (terminal_only == 0) % consider running cost as well
        for t = 1:T
            cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(t,range)' - h_d(range,1))' * (h_traj(t,range)' - h_d(range,1));
        end
    end      
    cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(T+1,range)' - h_d(range,1))' * (h_traj(T+1,range)' - h_d(range,1));

    fprintf('Iteration number: %d, Cost: %f \n', iter, cost(iter,1));
    
    U = U_new;
end 

save('optimal_controls.mat', 'U_new');

figure()
plot(cost)
title('cost vs iterations');

figure()
plot(U_new)

%%

figure()
plot(x,h_traj(1,:));
hold on 
plot(x, h_d, '-r', 'LineWidth',5);
hold off
legend('actual','desired');
title('starting profile');

figure()
plot(x,h_traj(end,:));
hold on 
plot(x, h_d, '-r', 'LineWidth',5);
hold off
legend('actual','desired');
title('end profile');

figure()
s = surf(x,(1:1:T+1),h_traj);
s.EdgeAlpha = 0.25;
xlabel('space')
ylabel('time')
zlabel('velocity')

ylim1 = min(min(h_traj));
ylim2 = max(max(h_traj));

fprintf("\n");
animate = input('Would you like to watch an animation of the field? Type 1: Yes or 0: No and press Enter !');

if animate
    fig = figure();
    for i = 1:T+1
        i
        plot(x, h_traj(i,:));
        hold on 
        plot(x, h_d, '-r');
        hold off
        legend('actual','desired');
        ylim manual
        ylim([ylim1 ylim2])
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        set(gcf, 'Toolbar', 'none', 'Menu', 'none');
        pause(0.01);
        clf(fig);
    end
end

% 
%  
% figure()
% plot(U_new)
% 
% mean_h = mean(h_traj, 1);
% var_h = var(h_traj,1); 
% figure()
% plot(x,h_d, '-r');
% hold on; 
% plot(x, mean_h,'-b');
% hold on; 
% plot(x, mean_h + 2.*sqrt(var_h),'--k');
% hold on; 
% plot(x, mean_h - 2.*sqrt(var_h),'--k');