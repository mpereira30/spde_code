
function [t,ut] = pde_fd_white_noise_series(u0,T,a,N,J,bctype, sigma, epsilon)

    Dt = T/N; 
    t = [0:Dt:T]'; 
    h = a/J;
    
    % set matrix A according to boundary conditions
    e = ones(J+1,1); 
    A = spdiags([-e 2*e -e], -1:1, J+1, J+1);
    
    switch lower(bctype)
        case {'dirichlet','d'}
            ind = 2:J; 
            A = A(ind,ind);
        case {'periodic','p'}
            ind = 1:J; 
            A = A(ind,ind); 
            A(1,end) = -1; 
            A(end,1) = -1;
        case {'neumann','n'}
            ind = 1:J+1; 
            A(1,2) = -2; 
            A(end,end-1) = -2;
    end
    
    EE = speye(length(ind))+Dt*epsilon*A/h/h;
    ut = zeros(J+1,length(t)); % initialize vectors
    ut(:,1) = u0; 
    u_n = u0(ind); % set initial condition
    
    bj = get_onedD_bj_white_noise(Dt,J,a);
    
    for k = 1:N, % time loop
%         if mod(k,100)==0
%             k
%         end
        nn = randn(length(bj),1);
        X = bsxfun(@times,bj,nn);
        dW = dst(X);
        u_new = EE\(u_n + sigma*dW);
        ut(ind, k+1) = u_new;
        u_n = u_new;
    end
    
    if lower(bctype)=='p' | lower(bctype)=='periodic'
        ut(end,:) = ut(1,:); % correct for periodic case
    end
    
end

