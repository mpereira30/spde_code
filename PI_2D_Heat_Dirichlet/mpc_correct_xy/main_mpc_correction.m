clc
clear 
close all

global J;
global N;
global a;
global mu_x;
global mu_y;
global sig_xx;
global sig_yy;
global T;
global dt;
global sigma;
global h_d;
global rollouts;
global rho;
global scale_factor;
global epsilon;
global terminal_only;
global sqr_loss; 
global r_x1;
global r_x2;
global r_y1;
global r_y2;

num_contours = 50;
colormap_val = jet;

terminal_only = 0; % Set 1 for considering only terminal state cost
sqr_loss = 1; 

% Set parameters:
dt = 0.01;

Tsim = 1.0; % in seconds
Tsim_steps = round(Tsim/dt);

Tf = 0.05; % final time (in seconds)
T = round(Tf/dt); % Horizon length

a = 0.5; % rod length 
J = 64; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';
y = (0:z:a)';
iters = 5; % number of PI iterations
rollouts = 25; % number of rollouts for sampling 

rho = 1000; % temperature parameter
sigma = 1/(sqrt(rho)); % standard deviation for Q-Wiener noise
scale_factor = 100;
epsilon = 1; % thermal diffusivity of the rod

%  Set initial and desired temperature distributions of rod:
desired_temp = 1;
h_d = NaN(J-1, J-1);
% h0 = zeros(J-1, J-1);
init_sigma = 0.5;
h0 = init_sigma * rand(J-1, J-1);

% Each column represents a patch. The elements of each column, are the
% lines at x1, x2, y1 and y2 whose intersection gives the patch.
r_y1 = round([0.48*(J-1), 0.48*(J-1), 0.48*(J-1), 0.18*(J-1), 0.78*(J-1)]);
r_y2 = round([0.52*(J-1), 0.52*(J-1), 0.52*(J-1), 0.22*(J-1), 0.82*(J-1)]);
r_x1 = round([0.48*(J-1), 0.18*(J-1), 0.78*(J-1), 0.48*(J-1), 0.48*(J-1)]);
r_x2 = round([0.52*(J-1), 0.22*(J-1), 0.82*(J-1), 0.52*(J-1), 0.52*(J-1)]);

% y-axis is row-wise
% x-axis is column-wise
h_d(r_y1(1):r_y2(1), r_x1(1):r_x2(1)) = 0.5*desired_temp; % desired temperature of center
for i = 2:length(r_x1)
   h_d(r_y1(i):r_y2(i), r_x1(i):r_x2(i)) = desired_temp;
end

% Actuator locations:
mu_y = [0.5; 0.2; 0.5; 0.8; 0.5].*a;
mu_x = [0.2; 0.5; 0.5; 0.5; 0.8].*a; 

sig_val = (0.1*a) * ones(length(mu_x),1);
sig_xx = sig_val.^2;
sig_yy = sig_val.^2;

N = length(mu_x); % number of actuators

% Initialize control sequences:
U = randn(T,N);

precompute = 1;

% First compute matricies curly-M, M and the curly-v-tilde vector offline:
if (precompute)
    fprintf('Computing curly M, capital M and curly v tilde\n');
    curly_M = compute_curly_M();
    M = compute_capital_M();
    curly_v_tilde = compute_curly_v_tilde();
    save('precompute.mat','M','curly_M','curly_v_tilde');
else
    load('precompute.mat');
    fprintf('Completed!. Starting iterations ....\n');
end
    
cost = zeros(Tsim_steps,1);
U_applied = zeros(Tsim_steps,N);
h_actual = zeros(J-1,J-1,Tsim_steps+1);
h_actual(:,:,1) = h0;

% For par-for (make local copies of gloval variables):
J_ = J; a_ = a; T_ = T; dt_ = dt; sigma_ = sigma; epsilon_ = epsilon; 

for cur_step = 1:Tsim_steps
    num_contours = 50;
    tic
    fprintf('current timestep: %d of %d,',cur_step, Tsim_steps );
    
    h0_ = h0(:);
    
%     cost_mpc = zeros(iters,1);
    for iter = 1:iters

        % Generate rollouts and noise samples:
        h_samples = zeros(rollouts,T+1,(J-1),(J-1)); % trajectories of spatial points
        xi_samples = zeros(rollouts,T,(J-1),(J-1)); % corresponding noise trajectories 
        
        parfor r = 1:rollouts
            [h_traj, xi_traj] = generate_rollouts(h0_, U, curly_v_tilde, 0, J_, a_, T_, dt_, sigma_, epsilon_);
            h_samples(r,:,:,:) = h_traj;
            xi_samples(r,:,:,:) = xi_traj;
        end

        % Cost computation and control update:
        U_new = PI_control(h_samples, xi_samples, U, curly_M, M);
        
        U = U_new;
    end

    %% Noise-free rollout and noise-free cost computation:
%     [h_traj, xi_traj] = generate_rollouts(h0_, U_new, curly_v_tilde, 1, J_, a_, T_, dt_, sigma_, epsilon_);
     
%     % With disturbances in the actual dynamics:
    sigma_sys = 0.2;
    [h_traj, xi_traj] = generate_rollouts(h0_, U_new, curly_v_tilde, 0, J_, a_, T_, dt_, sigma_sys, epsilon_);    
        
    h0 = squeeze(h_traj(2,:,:));
    h_actual(:,:,cur_step+1) = h0;
        
    % Set the initial control sequence for the next time step:
    U_applied(cur_step, :) = U_new(1,:);
    U = [U_new(2:end,:); U_new(end,:)];
    
    if(sqr_loss==1)
        for i = 1:length(r_x1)
            cost(cur_step, 1) = cost(cur_step, 1) + scale_factor * sum(sum( (h0(r_x1(i):r_x2(i), r_y1(i):r_y2(i)) - h_d(r_x1(i):r_x2(i), r_y1(i):r_y2(i))).^2 )); 
        end
    else
        for i = 1:length(r_x1)
            cost(cur_step, 1) = cost(cur_step, 1) + scale_factor * sum(sum( abs(h0(r_x1(i):r_x2(i), r_y1(i):r_y2(i)) - h_d(r_x1(i):r_x2(i), r_y1(i):r_y2(i))) )); 
        end
    end
    fprintf('cost: %f\n', cost(cur_step,1));

    toc
    fprintf('\n');
end

cmin = min(min(min(h_actual)));
cmax = max(max(max(h_actual)));

figure(1)
Z = padarray(h_d, [1 1], 0, 'both');
contourf(x,y,Z, num_contours);
colormap(colormap_val);
colorbar;
caxis([cmin cmax]);
title('desired profile');

figure(2)
plot(cost)

figure(3)
Z = squeeze(h_actual(:,:,1));
Z = padarray(Z, [1 1], 0, 'both');
%surf(x,y,Z);
contourf(x,y,Z, num_contours);
colormap(colormap_val);
colorbar;
caxis([cmin cmax]);
title('start');

figure(4)
Z = squeeze(h_actual(:,:,round(0.5*(Tsim_steps+1))));
Z = padarray(Z, [1 1], 0, 'both');
%surf(x,y,Z);
contourf(x,y,Z, num_contours);
colormap(colormap_val);
colorbar;
caxis([cmin cmax]);
title('half-way');

figure(5)
Z = squeeze(h_actual(:,:,end));
Z = padarray(Z, [1 1], 0, 'both');
%surf(x,y,Z);
contourf(x,y,Z, num_contours);
colormap(colormap_val);
colorbar;
caxis([cmin cmax]);
title('end');

figure(6)
plot(U_applied);
legend('act_1', 'act_2', 'center', 'act_4', 'act_5');

% 
% fig = figure();
% for i = 1:Tsim_steps
%     Z = squeeze(h_actual(:,:,i));
%     Z = padarray(Z, [1 1], 0, 'both');
%     contourf(x,y,Z, num_contours);
%     colormap(colormap_val)    
%     colorbar;
%     caxis([cmin cmax]);
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     pause(0.001);
%     clf(fig);
% end
