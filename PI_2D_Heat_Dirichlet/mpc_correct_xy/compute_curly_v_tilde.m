function [curly_v_tilde] = compute_curly_v_tilde()

%{  
    What is curly_v_tilde?

    * The infinite-dimensional controller is parameterized as,
      U(t, y_k, x_k) = m(y_k, x_k)'*u(t) = v(t, k_1, k_2) {whose shapes are : (1xN) x (Nx1)}
      v(t, k_1, k_2) is curly_v which is used in generate_rollouts().
    
    * But, before that, we need to compute m(y_k, x_k) for every (x_k,y_k). 
      curly_v_tilde is m(y_k, x_k) computed of each actuator n and each (x_k, y_k).
      After the last 2 dimensions are vectorized, we get the shape of { N x (J-1)^2 }
 
    * U(t, y_k, x_k) is the effect of the infinite-dimensional control at a particular spatial point (x_k, y_k) even if we have only a
      handful of actuators (eg. 5 and not J-1 x J-1). curly_v_tilde gives the appropriate scaling of the N controls at a particular spatial
      position depending on its co-ordinates.
  
%}
    global N;
    global a;
    global mu_x;
    global mu_y;
    global sig_xx;
    global sig_yy;
    global J;

    temp = zeros(N, (J-1), (J-1)); % shape is [ (number of actuators) x (y_dim) x (x_dim) ]
    curly_v_tilde = zeros(N, (J-1)^2); % (N x y_dim * x_dim)
    for n = 1:N 
       for i = 1:(J-1) % y_dim (Dirichlet B.C.s)
          for j = 1:(J-1) % x_dim (Dirichlet B.C.s)
             y = i*a/J; % y_dim <=> rows
             x = j*a/J; % x_dim <=> columns
             temp(n,i,j) = exp( -0.5 * ((x - mu_x(n))^2)/sig_xx(n) ) * exp( -0.5 * ((y - mu_y(n))^2)/sig_yy(n) ); 
          end
       end
       temp_1 = squeeze(temp(n,:,:));
       temp_1 = (temp_1(:))'; % equivalent to vec operation in linear algebra
       curly_v_tilde(n,:) = temp_1;
    end
    
end