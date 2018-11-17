function [curly_v_tilde] = compute_curly_v_tilde()

    global J;
    global N;
    global a;
    global mu;
    global sig;

    % generate the curly_v_tilde matrix:
    curly_v_tilde = zeros(N,J-1);
    for j = 1:(J-1) % iterate through spatial points
        m = zeros(N,1);
        for n = 1:N % iterate through actuators
            m(n,1) = exp( -(1/(2*sig(n))) * ((j*a/J) - mu(n))^2 );
            % ( This gives the effect of the nth actuator on the jth
            % spatial position )
        end
        curly_v_tilde(:,j) = m; 
    end


end