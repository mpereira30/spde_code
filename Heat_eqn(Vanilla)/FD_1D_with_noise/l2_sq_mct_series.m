
function [out] = l2_sq_mct_series(T,a,N,J,M,epsilon,sigma,bctype)
    
    v = 0;
    u0 = zeros(J+1,1);
    
    parfor i=1:M,
        [t, ut] = pde_fd_white_noise_series(u0,T,a,N,J,bctype,sigma,epsilon);
        v = v + sum(ut(1:end-1,end).^2); % cumulative sum (for Dirichlet BCs)
    end;

    out= v*a/J/M; % return average
end