function m = matrix_exp(A)
    [V, D, W] = eig(A);
    
    m = W * exp(D) * inv(W);
end