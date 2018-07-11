function a = A(WG, z)
    a = [];
       
    betas = WG.betas(:, z*WG.steps+1);
    
    s = size(betas,1);
    for m = 1:s
        for n = 1:m
            a(n,m) = Anm(WG, z, n, m);
            a(m,n) = -real(a(n,m)) + (imag(a(n,m)))*1i;
        end
    end
    
end