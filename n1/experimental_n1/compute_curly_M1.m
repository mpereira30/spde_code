function [curly_M] = compute_curly_M1()
    
global J;
global N;
global a;
global mu;
global sig;
    
%% Using Numerical integration:

    curly_M_num_int = zeros(J-1,N); % each row of this matrix is m^_i

    for j = 1:J-1 % by Marcus
       for i = 1:N 

           cur_mu = mu(i);
           cur_sig = sig(i); 
           fun = @(x, cur_sig, cur_mu) exp( -(1/(2*cur_sig)) .* (x - cur_mu).^2 ) .* sin( (pi*j/a) .* x);
           curly_M_num_int(j,i) =  sqrt(2/a) * integral(@(x) fun(x, cur_sig, cur_mu),0,a);

       end
    end
    curly_M = curly_M_num_int;

end