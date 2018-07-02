function a = A(WG, z)
    a = [];
    betas = WG.betas(:, z*20);
    profiles = WG.mode_profiles(z);
    
    s = size(betas,2);
    for m = 1:s
        for n = 1:m
            a(n,m) = Anm(WG, profiles, betas, z, n, m);
            a(m,n) = -real(a(n,m)) + imag(a(n,m));
        end
    end
    
end