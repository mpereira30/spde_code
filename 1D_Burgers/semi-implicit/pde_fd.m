
function [ut] = pde_fd(u0, dt, h, a, N, J, method, nu, dbc_val, diff_scheme, add_noise, sigma)
    
    if method == 's' % semi-implicit numerical scheme

        % set matrix A according to boundary conditions
        e       = ones(J+1,1); 
        A       = spdiags([-e 2*e -e], -1:1, J+1, J+1);
        ind     = 2:J; % Incorporating Dirichlet B.C.s
        A       = A(ind,ind);
        EE      = speye(length(ind)) + dt * nu * A/h^2;    
        
        root_q  = ones(J-1,1); % eigen values for cylindrical Q-Wiener 
        bj      = root_q * sqrt(2*dt/a);                  
        
        ut        = zeros(J+1, N+1); % Container to store the time evolution of the field
        ut(1,:)   = dbc_val(1)*ones(1,N+1);
        ut(end,:) = dbc_val(2)*ones(1,N+1);
        ut(:,1)   = u0; % overwrite with the initial profile 
        u_n       = u0(ind); % set profile at nth timestep 

        dbc_vec        = zeros(J-1,1);
        dbc_vec(1,1)   = dt * nu * dbc_val(1) / h^2;
        dbc_vec(J-1,1) = dt * nu * dbc_val(2) / h^2;
        
        for k = 1:N % time loop
            if mod(k,100)==0
                k
            end

            % advection term (non-linearity)
            if diff_scheme == 'c' % Using central difference
                fn = u_n .* (0.5 * dt / h) .* ( [u_n(2:end,1);dbc_val(2)] - [dbc_val(1);u_n(1:end-1,1)] ); 
            elseif diff_scheme == 'b' % Using backward difference
                fn = u_n .* (dt / h) .* ( u_n - [dbc_val(1); u_n(1:end-1,1)] ); 
            end
            
            % Generate space-time white noise:
            if add_noise == 1
                nn = randn(length(bj),1);
                X = bsxfun(@times,bj,nn);
                dW = dst(X);
            else 
                dW = zeros(size(u_n));
            end
            
            u_new = EE\(u_n -fn + sigma*dW + dbc_vec);
            ut(ind, k+1) = u_new;
            u_n = u_new;
        end
        
    elseif method == 'e' % fully explicit numerical scheme
        
        %%%%%%%%%%%%%%%% CAUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THIS IS INCOMPLETE, HAVE TO ACCOUNT FOR NON-ZERO B.Cs IN LAPLACIAN 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ut = zeros(J+1, N+1); % This takes care of homogeneous Dirichlet B.C.s at indices 1 and (J+1)
        ut(1,:) = dbc_val(1)*ones(1,N+1);
        ut(end,:) = dbc_val(2)*ones(1,N+1);
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

