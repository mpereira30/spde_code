
function [h_traj, xi_traj] = generate_rollouts(h0, U, curly_v_tilde, noise_free, J, nu, a, T, dbc_val, dt, sigma)
    
    z = a/J;
    
    % set matrix A according to boundary conditions
    e       = ones(J+1,1); 
    A       = spdiags([-e 2*e -e], -1:1, J+1, J+1);
    ind     = 2:J; % Incorporating Dirichlet B.C.s
    A       = A(ind,ind);
    EE      = speye(length(ind)) + dt * nu * A/z^2;    

    root_q  = ones(J-1,1); % eigen values for cylindrical Q-Wiener 
    bj      = root_q .* sqrt(2*dt/a);                  

    ht        = zeros(J+1, T+1); % Container to store the time evolution of the field
    ht(1,:)   = dbc_val(1)*ones(1,T+1);
    ht(end,:) = dbc_val(2)*ones(1,T+1);
    ht(:,1)   = h0; % overwrite with the initial profile 
    h_n       = h0(ind); % set profile at nth timestep 

    dbc_vec        = zeros(J-1,1);
    dbc_vec(1,1)   = dt * nu * dbc_val(1) / z^2;
    dbc_vec(end,1) = dt * nu * dbc_val(2) / z^2;

    % generate the curly_v matrix:
    curly_v = (U * curly_v_tilde)'; %  U x V~ = (T, N) x (N, J-1) = (T,J-1)^T = (J-1,T)

    Xi = zeros(J-1,T);
    for t = 1:T % time loop

        % advection term (non-linearity)
        fn = h_n .* (0.5 * dt / z) .* ( [h_n(2:end,1);dbc_val(2)] - [dbc_val(1);h_n(1:end-1,1)] ); % using central difference

        % Generate space-time white noise:
        if noise_free == 0
            Xi(:,t) = randn(J-1,1);
            dW = dst( bsxfun(@times,bj,Xi(:,t)) );
        else 
            dW = zeros(size(h_n));
        end
        
        h_new = EE\(h_n -fn + sigma*dW + dbc_vec + dt*curly_v(:,t));
        ht(ind, t+1) = h_new;
        h_n = h_new;
    end
    
    h_traj = ht';
    xi_traj = Xi';

end

