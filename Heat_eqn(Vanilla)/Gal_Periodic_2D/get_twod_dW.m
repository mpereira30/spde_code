
function [dW1,dW2]=get_twod_dW(bj,kappa,M)

    J=size(bj);
    if(kappa==1) % get xi_j
        nnr=randn(J(1),J(2),M); 
        nnc=randn(J(1),J(2),M);
    else % sum over kappa steps
        nnr=squeeze(sum(randn(J(1),J(2),M,kappa),4));
        nnc=squeeze(sum(randn(J(1),J(2),M,kappa),4));
    end
    
    nn2 = nnr + sqrt(-1)*nnc; 
    tmphat = bsxfun(@times,bj,nn2);
    tmp = ifft2(tmphat); 
    dW1 = real(tmp); 
    dW2 = imag(tmp);
end