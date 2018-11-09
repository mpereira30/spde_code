
function [t,u,ut] = spde_twod_Gal(u0,T,a,N,J,alpha)
    
    M = 1; % number of samples
    kappa = 1; 
    sigma = 0.001; % Q-Wiener noise sigma
    
    dtref = T/N; 
    Dt = kappa * dtref; 
    t = [0:Dt:T]'; 
    ut = zeros(J(1)+1, J(2)+1, N);
        
    % Set Lin Operator
    lambdax = 2*pi*[0:J(1)/2 -J(1)/2+1:-1]'/a(1);
    lambday = 2*pi*[0:J(2)/2 -J(2)/2+1:-1]'/a(2);
    [lambdaxx lambdayy] = meshgrid(lambday,lambdax);
    A = ( lambdaxx.^2 + lambdayy.^2);
    MM = 1*A; % epsilon = 1
    EE = 1./(1+Dt*MM);
    
    bj = get_twod_bj(dtref, J, a, alpha); % get noise coeffs
    u = repmat( u0(1:J(1),1:J(2)), [1,1,M]); % initial condition, last parameter M =1
    uh = repmat(fft2(u0(1:J(1),1:J(2))),[1,1,M]);

    % initialize
    uh1 = zeros(J(1),J(2),M);
    ut = zeros(J(1)+1,J(2)+1,N);
    ut(:,:,1) = u0;
    
    for n=1:N/kappa, % time loop
        dW = get_twod_dW(bj, kappa, M); % kappa and M are 1 
        gudWh = fft2(sigma.*dW);
        uh_new = bsxfun(@times, EE, (uh+gudWh)); % update u (semi-implicit propagation step) 
        u = real(ifft2(uh_new)); 
        ut(1:J(1),1:J(2),n+1) = u(:,:,end);
        uh = uh_new;
    end
    
    % make periodic
    u(J(1)+1,:,:) = u(1,:,:);
    u(:,J(2)+1,:) = u(:,1,:); 
    ut(J(1)+1,:,:) = ut(1,:,:); 
    ut(:,J(2)+1,:) = ut(:,1,:);
end

