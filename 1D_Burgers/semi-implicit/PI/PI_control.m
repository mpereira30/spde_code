%     h_samples = zeros(rollouts,T+1,J+1); % trajectories of spatial points
%     xi_samples = zeros(rollouts,T,J-1); % corresponding noise trajectories 

function [U_new, avg_cost] = PI_control(h_samples, xi_samples, U, curly_M, M)

    global J N a mu sig T dt sigma h_d rollouts rho scale_factor range nu terminal_only gamma
    
    J_h = zeros(rollouts,1);
    for r = 1:rollouts
        if(terminal_only == 0) % consider running cost as well. 
            for t = 1:T
                J_h(r,1) = J_h(r,1) + scale_factor * (squeeze(h_samples(r,t,range)) - h_d(range,1))' * (squeeze(h_samples(r,t,range)) - h_d(range,1));
            end
        end        
        J_h(r,1) = J_h(r,1) + scale_factor * (squeeze(h_samples(r,end,range)) - h_d(range,1))' * (squeeze(h_samples(r,end,range)) - h_d(range,1));
    end
    
%     avg
    
    zeta_1 = zeros(rollouts,1);   
    for r = 1:rollouts
        for t = 1:T
           m_temp = zeros(1,N); 
           for j = 1:(J-1) 
                m_temp = m_temp + curly_M(j,:) * xi_samples(r,t,j);
           end
           zeta_1(r,1) = zeta_1(r,1) + sqrt(dt/rho) * U(t,:) * m_temp'; % TODO: sqrt(dt) or dt ??
        end
    end

    zeta_2_temp = 0;
    for t = 1:T % goes to T only as controls are applied from timestep 0 to T-1 per algorithm 
        zeta_2_temp = zeta_2_temp + 0.5 * dt * U(t,:) * M * U(t,:)';
    end    
    zeta_2 = ones(rollouts,1) .* zeta_2_temp;
    
    % Compute total cost (with importance sampling term):
    J_h_tilde = J_h + zeta_1 + zeta_2;
    
    % Normalize the cost:
    minCost = min(J_h_tilde);
    maxCost = max(J_h_tilde);    
    J_h_tilde = J_h_tilde - (minCost .* ones(rollouts,1));
    J_h_tilde = J_h_tilde ./ (maxCost - minCost);
    
    weights = exp(-rho .* J_h_tilde);
    weights = weights ./ mean(weights);

    % Compute delta_u:
    U_new = zeros(T,N); 
    
    for t = 1:T
        
        u_update = zeros(1,N);
        for r = 1:rollouts
            m_temp = zeros(1,N); 
            for j = 1:(J-1) 
                m_temp = m_temp + curly_M(j,:) * xi_samples(r,t,j);
            end
            u_update = u_update + weights(r,1) * sqrt(dt) .* m_temp;
        end
        
        u_update = u_update/rollouts;         
        U_new(t,:) = U(t,:) + gamma * ( 1 / (dt*sqrt(rho)) ) .* (M \ u_update')';
    end
    
end















