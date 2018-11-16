clear; 
clc; 
close all;

dtref = 0.01; 
kappa = 1; 
r = 1/2; 
J = 128; 
a = 1;

bj = get_onedD_bj(dtref,J,a,r);

% discretized space:
T = 1;
t = 0:dtref:T;
x = zeros(J, length(t));
W = zeros(J, length(t));

% initialize the intermediate x points:
for i = 1:J-1
   x(i+1,:) = i*a/J; % we want x(1)=x(J+1)=0, so only initialize x(2) to x(J)
end

figure(1)
temp_t = zeros(J,length(t)); % for plotting purposes
for i = 2:1:length(t)
    temp_t(:,i) = t(i)*ones(J,1);
    W(2:J,i) = W(2:J,i-1) + get_onedD_dW(bj,kappa,0,1);
end
surf(x,temp_t,W);
grid on;
ylabel('time');
xlabel('x');
zlabel('W');
title('without kappa')

figure(2)
W = zeros(J, length(t));
temp_t = zeros(J,length(t)); % for plotting purposes
for i = 2:1:length(t)
    temp_t(:,i) = t(i)*ones(J,1);
    W(2:J,i) = W(2:J,i-1) + get_onedD_dW(bj,i-1,0,1);
end
surf(x,temp_t,W);
grid on;
ylabel('time');
xlabel('x');
zlabel('W');
title('with kappa')