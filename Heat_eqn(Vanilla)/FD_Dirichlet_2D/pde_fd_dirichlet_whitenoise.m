
% only dealing with Dirichlet boundary conditions:

function [t,ut] = pde_fd_dirichlet_whitenoise(u0, T, a, N, J, sigma)

    Dt=T/N;
    t=[0:Dt:T]'; 
    h=a/J;
    
    e = ones((J-1)*(J-1), 1);
    A = spdiags([-e -e 4*e -e -e], [-(J-1), -1, 0, 1, J-1], (J-1)*(J-1), (J-1)*(J-1));
        
    for i = 1:(J-1)-1
       A( i*(J-1)+1, i*(J-1) ) = 0;
       A( i*(J-1), i*(J-1)+1 ) = 0;
    end
    
    EE = speye(length(e)) + Dt*A/h^2;
    ut = zeros((J-1)*(J-1),length(t)); 
    ut(:,1) = u0;
    u_n = u0; % set initial condition
    
    for k = 1:N, % time loop
        if mod(k,100)==0
            k
        end
        Wn = sqrt(Dt/h) * randn( (J-1) * (J-1), 1);
        u_new = EE\(u_n + sigma*Wn);
        ut(:, k+1) = u_new;
        u_n = u_new;
    end
  
end