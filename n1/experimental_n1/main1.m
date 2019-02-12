clc
clear 
close all

global J N a mu sig T dt sigma h_d rollouts rho scale_factor range epsilon terminal_only gamma

monte_carlo_iters   = 128;
terminal_only       = 0; % Set 1 for considering only terminal state cost 
gamma               = 1.0; % step size
fprintf('\nUsing step size of %1.3f\n', gamma);
policy_type         = "open_loop/";
exp_type            = "suppress/";

%-------------------------------- Set parameters --------------------------
if exp_type == "accel/" 
    Tf          = 1.5; 
    fprintf('\nTime Horizon: %1.3f\n', Tf);
elseif exp_type == "suppress/"
    Tf          = 5.0;
    fprintf('\nTime Horizon: %1.3f\n', Tf);
end

dt              = 0.01;
T               = round(Tf/dt); % Horizon length

a               = 5.0; % axon length 
J               = 128; % number of spatial points
z               = a/J; % spatial discretization
x               = (0:z:a)';
epsilon         = 1; % thermal diffusivity of the rod

iters           = 200; % number of PI iterations
rollouts        = 100; % number of rollouts for sampling 
rho             = 10000;
sigma           = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor    = 10000;
%--------------------------------------------------------------------------

%-------------------------------- Set actuators --------------------------
range           = round(0.7*J):round(0.99*J);
mu              = [0.2 0.3 0.4 0.5 0.6 0.7 0.8].*a; % mean location of actuator in space
sig_val         = (0.1*a) * ones(length(mu));
sig             = sig_val.^2; % standard deviation of actuation in space
N               = length(mu); % number of actuators
%--------------------------------------------------------------------------

%-----------------------Choose experiment type-----------------------------
common_h0 = 1./(1+exp(-(2-x)/sqrt(2))); % set initial condition
h_d = zeros(size(x)); 

% Set targets:
if exp_type == "accel/" 
    h_d(range,1) = 1; % force end of axon to reach potential faster
    fprintf('Task: Accelerating voltage potential through an axon')
elseif exp_type == "suppress/"
    h_d(range,1) = 0; % suppress axon potential at the end
    fprintf('\nTask: Suprressing voltage at end of an axon\n')
end
%--------------------------------------------------------------------------

% First compute matricies curly-M, M and the curly-v-tilde vector offline:
curly_M = compute_curly_M1();
M = compute_capital_M();
curly_v_tilde = compute_curly_v_tilde1();
f = @(h) h.*(1-h).*(h+0.5);

cost = zeros(iters,1);   
h0 = common_h0; % Reset the initial profile every monte-carlo iteration

U = zeros(T,N); % Initialize control sequence

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

    %% Test rollout and compute cost:
    [h_traj, xi_traj] = generate_rollouts1(h0, U_new, curly_v_tilde, 0, f);

    for t = 1:T
        cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(t,range)' - h_d(range,1))' * (h_traj(t,range)' - h_d(range,1));
    end
    cost(iter,1) = cost(iter,1) + scale_factor * (h_traj(T+1,range)' - h_d(range,1))' * (h_traj(T+1,range)' - h_d(range,1));

    fprintf('Iteration number: %d, Cost: %f \n', iter, cost(iter,1));

    U = U_new;
end

% Perform the Monte-Carlo iterations:
Uoptimal        = U_new;
mc_cost         = zeros(monte_carlo_iters,1);
allEndProfiles  = zeros(monte_carlo_iters, length(x));

for mc_iter = 1: monte_carlo_iters

    h0 = common_h0; % Reset the initial profile every monte-carlo iteration
    
    [h_traj, xi_traj] = generate_rollouts1(h0, Uoptimal, curly_v_tilde, 0, f);
    for t = 1:T
        mc_cost(mc_iter,1) = mc_cost(mc_iter,1) + scale_factor * (h_traj(t,range)' - h_d(range,1))' * (h_traj(t,range)' - h_d(range,1));
    end
    mc_cost(mc_iter,1) = mc_cost(mc_iter,1) + scale_factor * (h_traj(T+1,range)' - h_d(range,1))' * (h_traj(T+1,range)' - h_d(range,1));
    
    fprintf('MC iteration: %d, total cost = %f\n', mc_iter, mc_cost(mc_iter,1));
    allEndProfiles(mc_iter,:) = h_traj(end,:);
end

if exp_type == "accel/" 
    allEndProfilesAccel = allEndProfiles;
    save('accel_profiles.mat', 'allEndProfilesAccel');
elseif exp_type == "suppress/"
    allEndProfilesSuppress = allEndProfiles;
    save('suppress_profiles.mat', 'allEndProfilesSuppress');
end

%     filepath = "monte_carlo_data/" + policy_type + exp_type + num2str(mc_iter)+".mat";
%     fprintf('MC-Iteration complete! Saving data at: %s\n', filepath);   
%     endprofile = h_traj(end,:);
%     costprofile = cost;
%     save(filepath, 'endprofile', 'costprofile')    


% figure()
% plot(cost)
% 
% 
% %%
% 
% figure()
% plot(x,h_traj(1,:));
% title('start');
% 
% figure()
% plot(x, h_traj(round(0.25*T),:) );
% title('quarter-way');
% 
% figure()
% plot(x, h_traj(round(0.5*T),:) );
% title('mid-way');
% 
% figure()
% plot(x,h_traj(end,:));
% hold on
% plot(x, h_d, 'g')
% title('end');
% 
% fig = figure();
% for i = 1:T
%     plot(x,h_d,'-r');
%     hold on;
%     plot(x,h_traj(i,:),'-b');
%     ylim manual
%     ylim([-0.125 1.25])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.02);
%     clf(fig);
% end

% figure()
% surf(x,(1:1:T+1),h_traj)
% 
% figure()
% plot(U_new)
% 
% figure()
% contourf(h_traj)