clc
clear 
close all

load('openloop_data.mat')

a = 1; % rod length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';

% Set parameters:
dt = 0.01;

% Total simulation time and timesteps:
Tsim = 1.0; % in seconds
Tsim_steps = round(Tsim/dt);
desired_temp = 5;

fig = figure();
F(Tsim_steps) = struct('cdata',[],'colormap',[]);
for i = 1:Tsim_steps
    plot(x,h_d, '-r');
    hold on;
    plot(x,h_traj(i,:), '-b');
    hold on;
    plot(0.2,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    plot(0.5,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    plot(0.8,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    legend('desired temperature','actual termperature','Actuators','Location','northeastoutside')
    title('Open-loop for Heat SPDE with Stochastic disturbances')
    xlabel('spatial position along rod')
    ylabel('temperature')
    ylim manual
    ylim([0 desired_temp*1.25])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    drawnow
    F(i) = getframe(fig);
    clf(fig);
end

movie2avi(F,'myavifile.avi','Compression','None','fps',15,'quality',1)

