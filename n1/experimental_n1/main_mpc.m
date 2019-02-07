clc
clear 
close all

global J N a mu sig T dt sigma h_d rollouts rho scale_factor range epsilon terminal_only gamma

monte_carlo_iters   = 128;
terminal_only       = 0; % Set 1 for considering only terminal state cost
gamma               = 0.1; % step size
fprintf('\nUsing step size of %1.3f\n', gamma);

%-------------------------------- Set parameters --------------------------
dt              = 0.01;
Tsim            = 5.0; % Total simulation time in seconds
Tsim_steps      = round(Tsim/dt); % timesteps
Tf              = 0.1; % % MPC final time horizon 
T               = round(Tf/dt); % mpc timesteps

a               = 5.0; % axon length 
J               = 128; % number of spatial points
z               = a/J; % spatial discretization
x               = (0:z:a)';
epsilon         = 1; % viscosity coefficient 

range           = round(0.7*J):round(0.99*J);
mu              = [0.2 0.3 0.4 0.5 0.6 0.7 0.8].*a; % mean location of actuator in space
sig_val         = (0.1*a) * ones(length(mu));
sig             = sig_val.^2; % standard deviation of actuation in space
N               = length(mu); % number of actuators

iters           = 5; % number of PI iterations
rollouts        = 100; % number of rollouts for sampling 
rho             = 10000;
sigma           = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor    = 10000;
%--------------------------------------------------------------------------

%-----------------------Choose experiment type-----------------------------
common_h0 = 1./(1+exp(-(2-x)/sqrt(2))); % set initial condition
h_d = zeros(size(x)); 

% Set targets:
% h_d(range,1) = 1; % force end of axon to reach potential faster
% exp_type = "accel/";
% fprintf('Task: Accelerating voltage potential through an axon')

h_d(range,1) = 0; % suppress axon potential at the end
exp_type = "suppress/";
fprintf('\nTask: Suprressing voltage at end of an axon\n')
%--------------------------------------------------------------------------

%---------------------------Precompute matrices ---------------------------
curly_M = compute_curly_M1();
M = compute_capital_M();
curly_v_tilde = compute_curly_v_tilde1();
%--------------------------------------------------------------------------

f = @(h) h.*(1-h).*(h+0.5); % non-linear function in Nagumo equation

for mc_iter = 1:monte_carlo_iters

    fprintf('\nMonte Carlo Iteration: %d\n', mc_iter);    
    
    h0 = common_h0; % Reset the initial profile every monte-carlo iteration
    U = zeros(T,N); % Re-initialize control sequence
    cost = zeros(Tsim_steps,1);

    % Container to store the trajectory of the actual system:
    h_actual = zeros(Tsim_steps+1, J+1);
    h_actual(1,:) = h0; % store the initial temperature profile

    for cur_step = 1:Tsim_steps

        fprintf('time: %1.2f s of %1.2f s, ',cur_step*dt, Tsim_steps*dt );

        % MPC iterations:
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
        [h_traj, xi_traj] = generate_rollouts1(h0, U_new, curly_v_tilde, 0, f);

        % Take the propagated state after application of 1st control of U_new:
        h0 = ( h_traj(2,:) )';
        h_actual(cur_step+1, :) = h0';

        % Set the initial control sequence for the next time step:
        U = [U_new(2:end,:); U_new(end,:)];

        cost(cur_step) = scale_factor * (h0(range,1) - h_d(range,1))' * (h0(range,1) - h_d(range,1));
        fprintf('cost: %f\n', cost(cur_step));
    end
    
    filepath = "monte_carlo_data/mpc/"+exp_type+num2str(mc_iter)+".mat";
    fprintf('Iteration complete! Saving data at: %s\n', filepath);   
    endprofile = h_actual(end,:);
    costprofile = cost;
    save(filepath, 'endprofile', 'costprofile')

end

% figure()
% plot(cost)
% 
% 
% %%
% 
% figure()
% plot(x,h_actual(1,:));
% title('start');
% 
% figure()
% plot(x, h_actual(round(0.25*Tsim_steps),:) );
% title('quarter-way');
% 
% figure()
% plot(x, h_actual(round(0.5*Tsim_steps),:) );
% title('mid-way');
% 
% figure()
% plot(x,h_actual(end,:));
% hold on
% plot(x, h_d, 'g')
% title('end');
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
% 
% figure()
% surf((1:1:Tsim_steps+1),x,h_actual')
% 
% figure()
% plot(U_new)
% 
% figure()
% contourf(h_actual)