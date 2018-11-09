
% only dealing with Dirichlet boundary conditions:

function [ut, vt] = pde_fd_dirichlet(u0, v0, h, N, J, dt, nu, dbcvalue, num_diff_scheme)
    
    if num_diff_scheme == 'c'
        fprintf("Using central difference for advection terms");
    elseif num_diff_scheme == 'b'
        fprintf("Using backward difference for advection terms");
    end
    
    ut = zeros((J-1)*(J-1), N+1); 
    vt = zeros((J-1)*(J-1), N+1);
    ut(:,1) = u0;
    vt(:,1) = v0;
    u_n = u0;
    v_n = v0;
    
    
    e = ones((J-1)*(J-1), 1);
    A = spdiags([-e -e 4*e -e -e], [-(J-1), -1, 0, 1, J-1], (J-1)*(J-1), (J-1)*(J-1));
        
    for i = 1:(J-1)-1
       A( i*(J-1)+1, i*(J-1) ) = 0;
       A( i*(J-1), i*(J-1)+1 ) = 0;
    end
    
    EE = speye(length(e)) + dt * nu * A/h^2; % h = delta_x = delta_y

    for k = 1:N, % time loop
        if mod(k,100)==0
            k
        end
        
        % To compute the convection terms, first convert back to 2D matrix:
        un_mat = ( vec2mat(u_n, J-1) )';
        vn_mat = ( vec2mat(v_n, J-1) )';        
        
        % Below un_mat and vn_mat do not include the boundaries:
        % change in x : column 
        % change in y : row 
        if num_diff_scheme == 'b'
            
            dudx = un_mat - [dbcvalue*ones(J-1,1), un_mat(:,1:end-1)]; % prepend a column 
            dudy = un_mat - [dbcvalue*ones(1,J-1); un_mat(1:end-1,:)]; % prepend a row
            dvdx = vn_mat - [dbcvalue*ones(J-1,1), vn_mat(:,1:end-1)];
            dvdy = vn_mat - [dbcvalue*ones(1,J-1); vn_mat(1:end-1,:)];
            
            u_new = EE\( u_n - (dt/h).*u_n.*dudx(:) - (dt/h).*v_n.*dudy(:) );
            v_new = EE\( v_n - (dt/h).*u_n.*dvdx(:) - (dt/h).*v_n.*dvdy(:) );            
        elseif num_diff_scheme == 'c'
            
            dudx = [un_mat(:, 2:end), dbcvalue*ones(J-1,1)] - [dbcvalue*ones(J-1,1), un_mat(:,1:end-1)];
            dudy = [un_mat(2:end, :); dbcvalue*ones(1,J-1)] - [dbcvalue*ones(1,J-1); un_mat(1:end-1,:)];
            dvdx = [vn_mat(:, 2:end), dbcvalue*ones(J-1,1)] - [dbcvalue*ones(J-1,1), vn_mat(:,1:end-1)];
            dvdy = [vn_mat(2:end, :); dbcvalue*ones(1,J-1)] - [dbcvalue*ones(1,J-1); vn_mat(1:end-1,:)];
            
            u_new = EE\( u_n - (0.5*dt/h).*u_n.*dudx(:) - (0.5*dt/h).*v_n.*dudy(:) );
            v_new = EE\( v_n - (0.5*dt/h).*u_n.*dvdx(:) - (0.5*dt/h).*v_n.*dvdy(:) );
        end
        
        ut(:, k+1) = u_new;
        vt(:, k+1) = v_new;        
        u_n = u_new;
        v_n = v_new;
    end
  
end