% Returns one element (row n, column m) of the matrix A
function a = Anm (WG, z, n, m)
    if m == n
        a = 0;
        return;
    end
   
    numerator = WG.dedz(z) * WG.k * overlap_integral(WG, z, n, m);
    denominator = WG.betas(m, z*WG.steps+1) - WG.betas(n, z*WG.steps+1);

    p = phase(WG, z, n, m);
    
    a = (numerator / denominator) * p;    
end