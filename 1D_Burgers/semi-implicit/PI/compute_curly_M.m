function [curly_M] = compute_curly_M()
    
global J;
global N;
global a;
global mu;
global sig;

 %% Using DST:
% 
%     curly_M_dst = zeros(J-1,N); % each row of this matrix is m^_i
%     x_j = zeros(J-1,1);
% 
%     for j=1:(J-1)
%        x_j(j,1) = j*a/J; 
%     end
%     
%     for i=1:N
%        x_n = zeros(J-1,1);
%        for j = 1:(J-1) 
%             x_n(j,1) = exp( -(1/(2*sig(i))) * (x_j(j,1) - mu(i))^2 );       
%        end
% %        curly_M_dst(:,i) = idst(x_n) .* sqrt(a*0.5);
%        curly_M_dst(:,i) = dst(x_n) .* sqrt(a*2)/J; % delta t is not considered in computation of inner products
%     end
%     curly_M = curly_M_dst;
    
    
%% Using Numerical integration:

    curly_M_num_int = zeros(J-1,N); % each row of this matrix is m^_i

    for j = 1:(J-1)
       for i = 1:N 

           cur_mu = mu(i);
           cur_sig = sig(i); 
           fun = @(x, cur_sig, cur_mu) exp( -(1/(2*cur_sig)) .* (x - cur_mu).^2 ) .* sin( (pi*j/a) .* x);
           curly_M_num_int(j,i) =  sqrt(2/a) * integral(@(x) fun(x, cur_sig, cur_mu),0,a);

       end
    end
    curly_M = curly_M_num_int;

end