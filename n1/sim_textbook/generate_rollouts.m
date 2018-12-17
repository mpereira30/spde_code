
function [h_traj, xi_traj] = generate_rollouts(h0, noise_free, f)

    global J;
    global a;
    global T;
    global dt;
    global sigma;
    global epsilon;
    
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

    root_q = ones(J-1,1); 
    bj = root_q .* sqrt(2*dt/a); % root_q = 1 for cylindrical noise     

    Xi = zeros(J+1,T);
        
    for t = 1:T % time loop -> (1,T) as per MATLAB indexing. But, this is actually (0, T-1) as per algorithm.  
        Xi(:,t) = randn(J+1,1);
        dW = dst( bsxfun(@times,bj,Xi(2:J,t)) );
        dW = [0;dW;0];
        
        if(noise_free == 1)
            h_new = EE\(h_n + dt * f(h_n));
        else
            h_new = EE\(h_n + dt * f(h_n) + sigma*dW);
        end
        
        ht(ind, t+1) = h_new;
        h_n = h_new;
    end
    
    h_traj = ht';
    xi_traj = Xi';
    
end

