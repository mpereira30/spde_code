
function [h_traj, xi_traj] = generate_rollouts(h0, U, curly_v_tilde, noise_free,  J, a, T, dt, sigma, epsilon)

    z = a/J;
    
    % set matrix A according to boundary conditions
    e = ones((J-1)*(J-1), 1);
    A = spdiags([-e -e 4*e -e -e], [-(J-1), -1, 0, 1, J-1], (J-1)^2, (J-1)^2);
    
    % Account for Dirichlet B.Cs:
    for i = 1:(J-1)-1
       A( i*(J-1)+1, i*(J-1) ) = 0;
       A( i*(J-1), i*(J-1)+1 ) = 0;
    end    
    
    EE = speye(length(e)) + dt*epsilon*A/z/z;
    h_n = h0; % shape is ( (J-1)*(J-1), 1 )
    
    % generate the curly_v matrix:
    curly_v = (U * curly_v_tilde)'; % (T,N) * (N,(J-1)^2) = (T, (J-1)*(J-1))' = ((J-1)*(J-1),T)

    b = 2*sqrt(dt)/a;
    Xi = zeros(T,(J-1),(J-1));
    h_traj = zeros(T+1, J-1, J-1);
    h_traj(1,:,:) = ( vec2mat(h0,J-1) )';
        
    for t = 1:T % time loop -> (1,T) as per MATLAB indexing. But, this is actually (0, T-1) as per algorithm.  
        Xi(t,:,:) = randn((J-1),(J-1)); % j2 x j1
        
        % Generate 2D noise using DST:
%         dW = zeros((J-1),(J-1)); % shape is y_dim * x_dim
%         for k2 = 1:(J-1) % going along the row dims (y-dim)
%            for j2 = 1:(J-1)
%                x_n = sin(pi*j2*k2/J) * reshape(Xi(t,j2,:), [1,J-1]);
%                dW(k2,:) = dW(k2,:) +  dst(x_n);
%            end
%         end
%         dW = b.*dW;

        dW = b .* dst(dst(squeeze(Xi(t,:,:)))')';
        
        if(noise_free==1)
            h_new = EE\(h_n + dt*curly_v(:,t));
        else
            h_new = EE\(h_n + dt*curly_v(:,t) + sigma*dW(:));
        end
        
        h_traj(t+1,:,:) = ( vec2mat(h_new,J-1) )';
        h_n = h_new;
    end
    
    xi_traj = Xi;
    
end

