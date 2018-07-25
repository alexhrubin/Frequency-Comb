function a = A(WG, betas, i)
    a = [];
       
    zq = betas(end, 1:i);
    
    s = size(betas,1) - 1;
    for m = 1:s
        for n = 1:m
            if (n == m)
                a(n,m) = 0;
                continue
            end
            
            beta_diff = betas(m, 1:i) - betas(n, 1:i);
            integral = trapz(zq, beta_diff);
            phase = exp(1j * integral);
            
            numerator = WG.dedz() * WG.k * overlap_numeric(WG, betas, i, m, n);
            denominator = betas(m, i) - betas(n, i);
            
            a(n,m) = phase * numerator / denominator;
            a(m,n) = -real(a(n,m)) + (imag(a(n,m)))*1i;
        end
    end
end