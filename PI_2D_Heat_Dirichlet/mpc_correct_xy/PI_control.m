
function U_new = PI_control(h_samples, xi_samples, U, curly_M, M)

global J;
global N;
global a;
global mu_x;
global mu_y;
global sig_xx;
global sig_yy;
global T;
global dt;
global sigma;
global h_d;
global rollouts;
global rho;
global scale_factor;
global epsilon;
global terminal_only;
global sqr_loss; 
global r_x1;
global r_x2;
global r_y1;
global r_y2;

    J_h = zeros(rollouts,1);
    for r = 1:rollouts
        
        % Compute running cost:
        if(terminal_only==0 && sqr_loss==1 ) % consider running cost as well. 
            for t = 1:T
                for i = 1:length(r_x1)
                    range_y = r_y1(i):r_y2(i);
                    range_x = r_x1(i):r_x2(i);
                    J_h(r,1) = J_h(r,1) +  scale_factor * sum(sum( (squeeze( h_samples(r,t,range_y,range_x) ) - h_d(range_y,range_x)).^2 ));
                end
            end
        elseif(terminal_only==0 && sqr_loss==0 )
            for t = 1:T
                for i = 1:length(r_x1)
                    range_y = r_y1(i):r_y2(i);
                    range_x = r_x1(i):r_x2(i);
                    J_h(r,1) = J_h(r,1) +  scale_factor * sum(sum( abs(squeeze( h_samples(r,t,range_y,range_x) ) - h_d(range_y,range_x)) ));
                end
            end           
        end
        
        % Compute terminal cost:
        if(sqr_loss==1)
            for i = 1:length(r_x1)
                range_y = r_y1(i):r_y2(i);
                range_x = r_x1(i):r_x2(i);
                J_h(r,1) = J_h(r,1) +  scale_factor * sum(sum( (squeeze( h_samples(r,T+1,range_y,range_x) ) - h_d(range_y,range_x)).^2 ));
            end    
        else
            for i = 1:length(r_x1)
                range_y = r_y1(i):r_y2(i);
                range_x = r_x1(i):r_x2(i);
                J_h(r,1) = J_h(r,1) +  scale_factor * sum(sum( abs(squeeze( h_samples(r,T+1,range_y,range_x) ) - h_d(range_y,range_x)) ));
            end  
        end
            
    end
    
    zeta_1 = zeros(rollouts,1);   
    for r = 1:rollouts
        for t = 1:T
           % curly_M : (J-1, J-1, N), xi_samples : (rollouts,T,(J-1),(J-1))
           m_temp = reshape( sum( sum( bsxfun(@times, curly_M, squeeze(xi_samples(r,t,:,:))), 2 ), 1), [1,N]); 
           zeta_1(r,1) = zeta_1(r,1) + sqrt(dt/rho) * U(t,:) * m_temp';
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
            m_temp = reshape( sum( sum( bsxfun(@times, curly_M, squeeze(xi_samples(r,t,:,:))), 2 ), 1), [1,N]); 
            u_update = u_update + weights(r,1) * sqrt(dt) .* m_temp;
        end
        
        u_update = u_update/rollouts;         
        U_new(t,:) = U(t,:) + ( 1 / (dt*sqrt(rho)) ) .* (M \ u_update')';
    end
    
end















