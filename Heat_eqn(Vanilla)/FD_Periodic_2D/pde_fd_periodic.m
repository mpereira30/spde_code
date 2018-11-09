
function [t,ut] = pde_fd_periodic(u0, T, a, N, J)

    Dt=T/N
    t=[0:Dt:T]'; 
    h=a/J;
    
    % For periodic:
    e = ones(J*J, 1);
    A = spdiags([-e -e -e 4*e -e -e -e], [-(J-1)*J ,-J, -1, 0, 1, J, (J-1)*J], J*J, J*J);

    for i = 1:J-1
       A( (i*J)+1, i*J ) = 0;
       A( i*J, (i*J)+1 ) = 0;

       A( i*J, (i-1)*J + 1 ) = -1;
       A( (i-1)*J + 1, i*J ) = -1;

    end
    A( J*J, (J-1)*J + 1 ) = -1;
    A( (J-1)*J + 1, J*J ) = -1;

    EE = speye(length(e)) + Dt*A/h^2;
    ut = zeros(J*J, length(t)); 
    ut(:,1) = u0;
    u_n = u0; % set initial condition
    
    for k = 1:N, % time loop
        if mod(k,100)==0
            k
        end
        u_new = EE\u_n;
        ut(:, k+1) = u_new;
        u_n = u_new;
    end
  
end