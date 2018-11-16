function dW=get_onedD_dW(bj,kappa,iFspace,M)

if(kappa==1) % generate xi_j
    nn=randn(length(bj),M);
else % sum over kappa steps
    nn=squeeze(sum(randn(length(bj),M,kappa),3));
end

X=bsxfun(@times,bj,nn);

if(iFspace==1) % return b_j xi_j
    dW=X;
else % return Wiener increments
    dW=dst(X);
    
end

