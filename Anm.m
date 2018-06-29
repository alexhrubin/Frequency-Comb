% Returns one element (row n, column m) of the matrix A
function a = Anm (WG, profiles, betas, z, n, m)
    if m == n
        a = 0;
        return;
    end
   
    numerator = WG.dedz(z) * WG.k * overlap_integral(profiles, n, m);
    denominator = betas(m) - betas(n);
  
    p = phase(WG, m, n, z);
    
    a = (numerator / denominator) * p;
end