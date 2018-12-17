function [curly_v_tilde] = compute_curly_v_tilde1()

    global J;
    global N;
    global a;
    global mu;
    global sig;

    % generate the curly_v_tilde matrix for Neumann:
    curly_v_tilde = zeros(N,J+1);
    for j = 0:J % iterate through spatial points
        m = zeros(N,1);
        for n = 1:N % iterate through actuators
            m(n,1) = exp( -(1/(2*sig(n))) * ((j*a/J) - mu(n))^2 );
            % ( This gives the effect of the nth actuator on the jth
            % spatial position )
        end
        curly_v_tilde(:,j + 1) = m; 
    end


end