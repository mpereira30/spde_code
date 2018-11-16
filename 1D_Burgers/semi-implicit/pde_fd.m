
function [ut] = pde_fd(u0, dt, h, a, N, J, method, nu, bctype, dbc_val, diff_scheme, add_noise, sigma)
    
    if method == 's' % semi-implicit numerical scheme

        % set matrix A according to boundary conditions
        e = ones(J+1,1); 
        A = spdiags([-e 2*e -e], -1:1, J+1, J+1);

        switch lower(bctype)
            case {'dirichlet','d'}
                ind = 2:J; 
                A = A(ind,ind);
                root_q = ones(J-1,1); 
                bj = root_q * sqrt(2*dt/a); % root_q = 1 for cylindrical noise                 
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
            
            % Generate space-time white noise:
            if add_noise == 1
                nn = randn(length(bj),1);
                X = bsxfun(@times,bj,nn);
                dW = dst(X);
            else 
                dW = zeros(size(u_n));
            end
            
            u_new = EE\(u_n -fn + sigma*dW);
            ut(ind, k+1) = u_new;
            u_n = u_new;
        end

        if lower(bctype)=='p' | lower(bctype)=='periodic'
            ut(end,:) = ut(1,:); % correct for periodic case
        end
        
    elseif method == 'e' % fully explicit numerical scheme
        
        ut = dbc_val * ones(J+1, N+1); % This takes care of homogeneous Dirichlet B.C.s at indices 1 and (J+1)
        ind = 2:J;
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
            
            lapacian =  (nu * dt / h^2) * ([u_n(2:end,1);dbc_val] - 2 .* u_n + [dbc_val;u_n(1:end-1,1)]);
            
            % Generate space-time white noise:
            if add_noise == 1
                nn = randn(length(bj),1);
                X = bsxfun(@times,bj,nn);
                dW = dst(X);
            else 
                dW = zeros(size(u_n));
            end
            
            u_new = u_n -fn + lapacian + sigma*dW;
            ut(ind, k+1) = u_new;
            u_n = u_new;
        end      
        
    end
    
end

