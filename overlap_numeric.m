function overlap = overlap_numeric(WG, z, m, n, resolution)

    if (mod(m,2)==1 && mod(n,2)==0)
        overlap = 0;
        return
    elseif (mod(m,2)==0 && mod(n,2)==1)
        overlap = 0;
        return
    end
    
   mode_m = WG.eigenmode_function(z, m);
   mode_n = WG.eigenmode_function(z, n);
   product = @(x) mode_m(x) * mode_n(x);
   d = WG.width(z);
   X = -10*d:resolution:10*d;
   y = arrayfun(product, X);
   overlap = trapz(X, y);
   
end 