function [curly_M] = compute_curly_M1()
    
global J;
global N;
global a;
global mu;
global sig;
    
%% Using Numerical integration:

    curly_M_num_int = zeros(J+1,N); % each row of this matrix is m^_i
    
    % The index j stands for number of terms to approximate the Q-Wiener
    % increments. In the equation (10.16), they consider (J-1) terms to 
    % both discretize the spatial domain and (J-1) eigenfunctions
    % to approximate the Q-Wiener increments. However, these need not be
    % the same number. Below, I have used (J+1) terms to approximate the
    % Q-Wiener increments. The index j cannot be zero as sin(0)=0. So this
    % index always begins with 1, but can go upto (J-1) or any other
    % number. However, whatever the number here has to be the same number
    % of noise values sampled. 
        
    for j = 1:J+1 % by Marcus
       for i = 1:N 

           cur_mu = mu(i);
           cur_sig = sig(i); 
           fun = @(x, cur_sig, cur_mu) exp( -(1/(2*cur_sig)) .* (x - cur_mu).^2 ) .* sin( (pi*j/a) .* x);
           curly_M_num_int(j+1,i) =  sqrt(2/a) * integral(@(x) fun(x, cur_sig, cur_mu),0,a);

       end
    end
    curly_M = curly_M_num_int;

end