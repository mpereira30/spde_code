clc
clear
close all

x = (-5:0.1:5);

mu = 0;
sigma = 2;

fx_n = zeros(length(x),1);
fx = zeros(length(x),1);

for i = 1:length(x)
    fx(i) = exp( -(x(i)^2)/(2*sigma^2) );
    fx_n(i) = fx(i) / (sigma*sqrt(2*pi));
    
end

figure()
plot(x,fx)
title('Unnormalized')

figure()
plot(x,fx_n)
title('normalized')