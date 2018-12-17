
function [h_traj, xi_traj] = generate_rollouts1(h0, U, curly_v_tilde, noise_free, f)

    global J N a mu sig T dt sigma epsilon
    
    z = a/J;
    
    % set matrix A according to boundary conditions
    e = ones(J+1,1); 
    A = spdiags([-e 2*e -e], -1:1, J+1, J+1);
    ind=1:J+1; 
    A(1,2)=-2; 
    A(end,end-1)=-2;
    
    EE = speye(length(ind)) + dt*epsilon*A/z/z;
    ht = zeros(J+1,T+1); % initialize vectors
    ht(:,1) = h0; 
    h_n = h0(ind); % set initial condition
    
    % generate the curly_v matrix:
    curly_v = zeros(J+1,T);
    for t = 1:T 
        curly_v(:,t) = (U(t,:) * curly_v_tilde)'; %  U_t x V~ = (1, N) x (N, J-1)
    end

    root_q = ones(J-1,1); % for Dirichlet boundary conditions only
    bj = root_q .* sqrt(2*dt/a); % root_q = 1 for cylindrical noise     

    Xi = zeros(J+1,T);
        
    for t = 1:T % time loop -> (1,T) as per MATLAB indexing. But, this is actually (0, T-1) as per algorithm.  
        Xi(:,t) = randn(J+1,1); % This (J+1) has nothing to do with number of spatial points or this being Neumann problem. This is simply the number of 
        % terms (noise and eigen functions) chosen to approximate the
        % Q-Wiener increments. 
        
        dW = dst( bsxfun(@times,bj,Xi(2:J,t)) );
        dW = [0;dW;0]; 
        
        % Why are the dW's zero at the terminal points in Neumann B.Cs?
        % Answer:
        % They are always zero, irrespective of the B.Cs. To see why this
        % is the case, refer to equation (10.16) of the computational SPDE
        % textbook. For the domain D=(0,a), we sample points to evaluate
        % Xhi's (eigenfunctions) for the dWs. We only sample the x_k's for k = 1, 2, ..., J-1
        % The reason we don't consider k = 0 and k = J, is because Xhi is
        % zero for these values. Therefore, for k = 0 and k = J, dW = 0. 
        
        % For the Dirichlet case, we know what the value of the field is at the boundary
        % points and so we don't have to solve for it at every timestep. In
        % case of Neumann B.Cs, we have to solve for the field at the
        % boundary points at every timestep. Therefore, we need to consider
        % what the dWs are at the boundary points for the Neumann B.Cs
        % case. However, the dW's are zeros because of the eigenfunctions.
        
        if(noise_free==1)
            h_new = EE\(h_n + dt*curly_v(:,t) + dt * f(h_n));
        else
            h_new = EE\(h_n + dt*curly_v(:,t) + dt * f(h_n) + sigma*dW);
        end
        
        ht(ind, t+1) = h_new;
        h_n = h_new;
    end
    
    h_traj = ht';
    xi_traj = Xi';
    
end

