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
Tsim = 1.0; % in seconds
Tsim_steps = round(Tsim/dt);

% MPC final time horizon and timesteps
Tf = 0.25; % in seconds
T = round(Tf/dt); 

a = 1; % rod length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
iters = 50; % number of PI iterations per MPC timestep
rollouts = 200; % number of rollouts for sampling 

rho = 1; % temperature parameter
sigma = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor = 100;
epsilon = 1; % thermal diffusivity of the rod

%  Set initial and desired temperature distributions of rod:
init_temp = 0;
desired_temp = 5;

h_d = zeros(J+1,1);

% 1.) Heat rod center
% range = round(0.45*J):round(0.55*J);
% h_d(range,1) = desired_temp;

% 2.) Heat 2 sides 
% range1 = round(0.15*J):round(0.25*J);
% range2 = round(0.75*J):round(0.85*J);
% range = [range1, range2];
% h_d(range1,1) = +desired_temp;
% h_d(range2,1) = +desired_temp;

% 3.) Heat 3 sides 
range1 = round(0.18*J):round(0.22*J);
range2 = round(0.48*J):round(0.52*J);
range3 = round(0.78*J):round(0.82*J);
range = [range1, range2, range3];
h_d(range1,1) = +desired_temp;
h_d(range2,1) = +desired_temp * 0.5;
h_d(range3,1) = +desired_temp;

% 4.) full range
% range = 1:(J+1);
% h_d(round(0.45*J):round(0.55*J),1) = desired_temp;

% 1.) Single acutator
% mu = [0.5]; % mean location of actuator in space
% sig_val = 0.2;
% sig = [sig_val].^2; % standard deviation of actuation in space

% 2.) Multiple acutators
% mu = [0.15, 0.18, 0.2, 0.22, 0.25, 0.3, 0.35, 0.4, 0.45, 0.48, 0.5, 0.52, 0.55, 0.6, 0.65, 0.7, 0.75, 0.78, 0.8, 0.82, 0.85]; % mean location of actuator in space
% mu = [0.2 0.3 0.4 0.5 0.6 0.7 0.8]; % mean location of actuator in space
mu = [0.2 0.5 0.8];
% mu = [0.2 0.8];
sig_val = 0.1 * ones(length(mu));
sig = sig_val.^2; % standard deviation of actuation in space

N = length(mu); % number of actuators

% Set the initial temperature profile of the rod:
% h0 = h_d;
h0 = zeros(J+1,1);

% First compute matricies curly-M, M and the curly-v-tilde vector offline:
curly_M = compute_curly_M();
M = compute_capital_M();
curly_v_tilde = compute_curly_v_tilde();

% Container to store the trajectory of the actual system:
h_actual = zeros(J+1, Tsim_steps+1);
h_actual(:,1) = h0; % store the initial temperature profile

% Initial control sequence:
U = randn(T,N);

for cur_step = 1:Tsim_steps
    
    fprintf('current timestep: %d of %d, ',cur_step, Tsim_steps );
    
    %% MPC iterations:
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
        U = U_new;
       
    end
    
    %% Play the control on system:
%     [h_traj, xi_traj] = generate_rollouts(h0, U_new, curly_v_tilde, 1);
    % shape of h_traj will be (T+1, J+1). This involves h0 as well. 
    
    % Play on NOISY DYNAMICS:
    [h_traj, xi_traj] = generate_rollouts(h0, U_new, curly_v_tilde, 0);
    
    % Take the propagated state after application of 1st control of U_new:
    h0 = ( h_traj(2,:) )';
    h_actual(:,cur_step+1) = h0;
    
    % Set the initial control sequence for the next time step:
    U = [U_new(2:end,:); U_new(end,:)];
    
    cost(cur_step) = scale_factor * (h0(range,1) - h_d(range,1))' * (h0(range,1) - h_d(range,1));
    fprintf('cost: %f\n', cost(cur_step));
        
end

figure()
plot(cost)

figure()
plot(x,h_actual(:,1));
title('start');

figure()
plot(x, h_actual(:,round(0.25*T)) );
title('quarter-way');

figure()
plot(x, h_actual(:,round(0.5*T)) );
title('mid-way');

figure()
plot(x,h_actual(:,end));
title('end');
% 
% fig = figure();
% for i = 1:Tsim_steps
%     plot(x,h_d, '-r');
%     hold on;
%     plot(x,h_actual(:,i), '-b');
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
surf(x,(1:1:Tsim_steps+1),h_actual')

 
figure()
plot(U_new)

mean_h = mean(h_actual, 2);
var_h = var(h_actual,0,2); 
figure()
plot(x,h_d, '-r');
hold on; 
plot(x, mean_h,'-b');
hold on; 
plot(x, mean_h + sqrt(var_h),'--k');
hold on; 
plot(x, mean_h - sqrt(var_h),'--k');