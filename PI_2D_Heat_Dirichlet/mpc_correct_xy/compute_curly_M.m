function [curly_M] = compute_curly_M()
    
global N;
global a;
global mu_x;
global mu_y;
global sig_xx;
global sig_yy;
global J;
    
    
%% Using Numerical integration:

    curly_M_num_int = zeros(J-1, J-1, N); 
    for j2 = 1:(J-1)
        for j1 = 1:(J-1)
            for n = 1:N
                  
                  % Compute as a double integral (NOT FEASIBLE- TAKES VERY LONG!!!): 
%                 fun = @(x,y) exp( -0.5 .* ((x - mu_x(n)).^2)./sig_xx(n) ) .* exp( -0.5 .* ((y - mu_y(n)).^2)./sig_yy(n) ) .* sin( (j1*pi/a) .* x ) .* sin( (j2*pi/a) .* y );
%                 curly_M_num_int(j1,j2,n) = (2/a) * integral2(fun, 0, a, 0, a);

                  % Compute as a product of 2 integrals:
                  fun_x = @(x) exp( -0.5 .* ((x - mu_x(n)).^2)./sig_xx(n) ) .* sin( (j1*pi/a) .* x ); 
                  fun_y = @(y) exp( -0.5 .* ((y - mu_y(n)).^2)./sig_yy(n) ) .* sin( (j2*pi/a) .* y ); 
                  curly_M_num_int(j2,j1,n) = (2/a) * integral(fun_x, 0, a) * integral(fun_y, 0, a);
            end
        end
    end
    curly_M = curly_M_num_int;

end