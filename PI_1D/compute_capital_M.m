
function [M] = compute_capital_M()

global J;
global N;
global a;
global mu;
global sig;

%% Using Numerical integration:

    M = zeros(N,N);
    
    for i = 1:N
       for j = 1:N
           
           mu_i = mu(i);
           mu_j = mu(j);
           sig_i = sig(i);
           sig_j = sig(j);
           fun = @(x, mu_i, mu_j, sig_i, sig_j) exp( -(1/(2*sig_i)) .* (x - mu_i).^2 ) .* exp( -(1/(2*sig_j)) .* (x - mu_j).^2 );
           M(i,j) = integral(@(x) fun(x, mu_i, mu_j, sig_i, sig_j), 0, a);
       end
    end
    
end



