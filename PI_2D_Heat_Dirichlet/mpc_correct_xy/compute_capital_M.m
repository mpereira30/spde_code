
function [M] = compute_capital_M()

global N;
global a;
global mu_x;
global mu_y;
global sig_xx;
global sig_yy;


%% Using Numerical integration:

    M = zeros(N,N);
    
    for i = 1:N
       for j = 1:N
           
           mu_x_i = mu_x(i);
           mu_y_i = mu_y(i);
           mu_x_j = mu_x(j);
           mu_y_j = mu_y(j);
            
           sig_xx_i = sig_xx(i);
           sig_yy_i = sig_yy(i);
           sig_xx_j = sig_xx(j);
           sig_yy_j = sig_yy(j);
           
           % In scalar notation (faster method):
           fun = @(x,y) exp( -0.5 .* ((x - mu_x_i).^2)./sig_xx_i ) .* exp( -0.5 .* ((y - mu_y_i).^2)./sig_yy_i ) .* exp( -0.5 .* ((x - mu_x_j).^2)./sig_xx_j ) .* exp( -0.5 .* ((y - mu_y_j).^2)./sig_yy_j );
           M(i,j) = integral2(fun, 0, a, 0, a);

           % Using product of 2 integrals:
%           fun_x = @(x) exp( -0.5 .* ((x - mu_x_i).^2)./sig_xx_i ) .* exp( -0.5 .* ((x - mu_x_j).^2)./sig_xx_j ); 
%           fun_y = @(y) exp( -0.5 .* ((y - mu_y_i).^2)./sig_yy_i ) .* exp( -0.5 .* ((y - mu_y_j).^2)./sig_yy_j ); 
%           M(i,j) = integral(fun_x, 0, a) * integral(fun_y, 0, a);           

           % In vector/matrix notation (Slower method):
%            fun = @(x, y) exp( -0.5 * ([x;y] - [mu_x_i;mu_y_i])' * inv([sig_xx_i, 0; 0, sig_yy_i]) * ([x;y] - [mu_x_i;mu_y_i]) ) * ...
%                          exp( -0.5 * ([x;y] - [mu_x_j;mu_y_j])' * inv([sig_xx_j, 0; 0, sig_yy_j]) * ([x;y] - [mu_x_j;mu_y_j]) );
%            
%            M(i,j) = integral2(@(x,y) arrayfun(fun, x, y), 0, a, 0, a);

       end
    end
    
end