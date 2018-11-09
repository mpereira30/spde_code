
function [out] = l2_sq_mct(T,a,N,J,M,sigma)
    
    v = 0;
    u0 = zeros((J-1)*(J-1),1);
    
    parfor i=1:M,
        [t, ut] = pde_fd_dirichlet_whitenoise(u0, T, a, N, J, sigma);
        v = v+sum(ut(:,end).^2); % cumulative sum
    end;

%     h_sq = (a/J)^2;
%     out= v*h_sq/M; % return average, equation is (1/M)*h*v_avg
    out= v*(a/J)/M;
end