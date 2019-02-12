clc 
close all 
clear 

Tf              = 1.0; % final time (in seconds)
dt              = 0.01;
T               = round(Tf/dt); % Horizon length or total number of timesteps
a               = 2; % rod length 
J               = 128; % number of spatial points
z               = a/J; % spatial discretization
x               = (0:z:a)';
num_trials      = 128;

desired_vel     = 2.0;
h_d             = NaN(J+1,1);
range1          = round(0.18*(J+1)):round(0.22*(J+1));
range2          = round(0.48*(J+1)):round(0.52*(J+1));
range3          = round(0.78*(J+1)):round(0.82*(J+1));
h_d(range1,1)   = desired_vel;
h_d(range2,1)   = desired_vel * 0.5;
h_d(range3,1)   = +desired_vel;

policy_type     = "mpc/";
% extract monte-carlo data:
all_end_profiles = zeros(num_trials, length(x));

filepath    = policy_type;
fprintf("\nextracting data from %s ...\n", filepath);

for trial = 1:num_trials
    load(filepath + int2str(trial) + ".mat");
    all_end_profiles(trial, :) = endprofile;
end

mean_mpc   = mean(all_end_profiles, 1);
std_mpc    = std(all_end_profiles,0,1);

ylim_lower = -0.0;
figure()
boundedline(x, mean_mpc', std_mpc', 'alpha','-r');
hold on;

load('../openloop_mc_profiles.mat');
boundedline(x, mean(allEndProfiles, 1)', std(allEndProfiles, 0, 1)', 'alpha','-b');
hold on;

plot(x, h_d, '-r', 'LineWidth',5);
hold on

sig = (0.1*a);

plot(0.2*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.2*a - sig, 0.2*a + sig], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;

plot(0.3*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.3*a - sig, 0.3*a + sig], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;

plot(0.5*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.5*a - sig, 0.5*a + sig], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;

plot(0.7*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.7*a - sig, 0.7*a + sig], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;

plot(0.8*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.8*a - sig, 0.8*a + sig], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');

legend('actual end profile','desired profile', 'actuator locations', 'actuator effect (1-sigma)','','');
title('end profile');
xlabel('spatial position in meters');
ylabel('velocity in m/s');
