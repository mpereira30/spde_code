% clc
% clear 
% close all

load('Accel_Nagumo.mat')

a = 10; % rod length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';

% Set parameters:
dt = 0.01;

% Total simulation time and timesteps:
Tsim = 2.5; % in seconds
Tsim_steps = round(Tsim/dt);
desired_temp = 5;
range = round(0.7*J):round(0.99*J);
h_d = zeros(length(x),1);
h_d(range,1) = 1;

fig = figure();
F(Tsim_steps) = struct('cdata',[],'colormap',[]);
for i = 1:Tsim_steps
    plot(x,h_d, '-r');
    hold on;
    plot(x,h_traj(i,:), '-b');
    hold on;
    plot(0.2*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    plot(0.3*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    plot(0.4*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;    
    plot(0.5*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    plot(0.6*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    plot(0.7*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;    
    plot(0.8*a,0,'-go','MarkerSize',15, 'MarkerFaceColor','g');
    hold on;
    legend('desired voltage','actual voltage','Actuators','Location','northeastoutside')
%     legend('actual voltage')
%     title('Propagation of voltage in an axon given by the deterministic Nagumo PDE')
    title('Open-loop for Nagumo SPDE for Acceleration of voltage propagation in axon ')
    xlabel('spatial position along axon')
    ylabel('voltage')
    ylim manual
    ylim([-0.125 1.25])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    drawnow
    F(i) = getframe(fig);
    clf(fig);
end

movie2avi(F,'myavifile.avi','Compression','None','fps',15,'quality',1)

