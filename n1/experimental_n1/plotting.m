clc
clear
close all

load('openloop_accel.mat')
load('openloop_supress.mat')

a = 10; % axon length 
J = 128; % number of spatial points
z = a/J; % spatial discretization
x = (0:z:a)';

% figure()
% plot(cost_accel)
% hold on 
% plot(cost_suppress)
% hold off
% legend('acceleration', 'suppression')

h_d_accel = NaN(size(x));
h_d_supp = NaN(size(x));
range = round(0.7*J):round(0.99*J);
h_d_accel(range,1) = 1; 
h_d_supp(range,1) = 0; 

ylim_lower = -0.1;

figure()
plot(x, final_accel_profile, 'r')
hold on 
plot(x, h_d_accel, '--r')
hold on
plot(x, suppress_end_profile, 'b')
hold on
plot(x, h_d_supp, '--b')
hold on;
plot(0.2*a, ylim_lower,'-go','MarkerSize',15);
hold on;
xtemp = [0.2*a - 0.1*a: 0.01: 0.2*a + 0.1*a];
norm = normpdf(x,0,1);
% line([0.2*a - 0.1*a, 0.2*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;
plot(0.3*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.3*a - 0.1*a, 0.3*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;
plot(0.4*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.4*a - 0.1*a, 0.4*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;    
plot(0.5*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.5*a - 0.1*a, 0.5*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;
plot(0.6*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.6*a - 0.1*a, 0.6*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;
plot(0.7*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.7*a - 0.1*a, 0.7*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;    
plot(0.8*a, ylim_lower,'-go','MarkerSize',15);
hold on;
line([0.8*a - 0.1*a, 0.8*a + 0.1*a], [ylim_lower, ylim_lower], 'Color','red', 'Marker', 'o');
hold on;
ylim([ylim_lower 1.25])
legend('acceleration actual', 'acceleration desired', 'suppression actual', 'suppression desired', 'actuator centers', 'actuator effect (1-sigma)')


