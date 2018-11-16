
function [ut] = pde_fd(u0, dt, h, N, J, method, nu, bctype, dbc_val, diff_scheme)

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
    
    EE = speye(length(ind)) + dt * nu * A/h^2;
    % Container to store the time evolution of the field:
    ut = dbc_val * ones(J+1, N+1); % This takes care of homogeneous Dirichlet B.C.s at indices 1 and (J+1)
    ut(:,1) = u0; % overwrite with the initial profile 
    
    u_n = u0(ind); % set profile at nth timestep 
    
    for k = 1:N % time loop
        if mod(k,100)==0
            k
        end
        
        % advection term (non-linearity)
        if diff_scheme == 'c' % Using central difference
            fn = u_n .* (0.5 * dt / h) .* ( [u_n(2:end,1);dbc_val] - [dbc_val;u_n(1:end-1,1)] ); 
        elseif diff_scheme == 'b' % Using backward difference
            fn = u_n .* (dt / h) .* ( u_n - [dbc_val; u_n(1:end-1,1)] ); 
        end
        
        u_new = EE\(u_n -fn);
        ut(ind, k+1) = u_new;
        u_n = u_new;
    end
    
    if lower(bctype)=='p' | lower(bctype)=='periodic'
        ut(end,:) = ut(1,:); % correct for periodic case
    end
    
end

