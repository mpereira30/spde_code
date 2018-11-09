
function bj = get_onedD_bj_white_noise(dtref, J, a)
    root_q = ones(J-1,1); % for Dirichlet boundary conditions 
    bj = root_q * sqrt(2*dtref/a); % root_q = 1 for cylindrical noise 
end
