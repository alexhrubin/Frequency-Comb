% This is the integral of dS E_m*E_n
function integral = overlap_integral(WG, betas, i, m, n)

    d = WG.width(betas(end, i));

    beta_m = betas(m, i);
    beta_n = betas(n, i);


    u_n = sqrt( (WG.n*WG.k)^2 - beta_n^2);
    w_n = sqrt( (beta_n^2 - WG.k^2));
    
    u_m = sqrt( (WG.n*WG.k)^2 - beta_m^2);
    w_m = sqrt( (beta_m^2 - WG.k^2));
    
        
    if (mod(n,2) == 0) && (mod(m,2) ==1) %symm * asymm
        integral = 0;
        return
    
    elseif (mod(m,2) == 0) && (mod(n,2) ==1) %symm * asymm
        integral = 0;
        return
     
        
    elseif (mod(n,2)==1) && (mod(m,2)==1) %symm * symm
        norm_n = 1/sqrt( (cos(d*u_n)^2)/w_n + d + sin(2*d*u_n)/(2*u_n));
        norm_m = 1/sqrt( (cos(d*u_m)^2)/w_m + d + sin(2*d*u_m)/(2*u_m));
        
        outside_wg = 2*cos(u_n*d)*cos(u_m*d) / (w_n+w_m);
        inside_wg = sin((u_m + u_n)*d) / (u_m+u_n) + sin((u_m-u_n)*d) / (u_m-u_n);
        
        integral = norm_n * norm_m * (outside_wg + inside_wg);
        return
        
 
    elseif (mod(n,2)==0) && (mod(m,2)==0) %asymm * asymm
        norm_n = 1/sqrt( (sin(d*u_n)^2)/w_n + d - sin(2*d*u_n)/(2*u_n));
        norm_m = 1/sqrt( (sin(d*u_m)^2)/w_m + d - sin(2*d*u_m)/(2*u_m));
        
        outside_wg = 2*sin(u_n*d)*sin(u_m*d) / (w_n+w_m);
        inside_wg = sin((u_m - u_n)*d) / (u_m-u_n) - sin((u_m+u_n)*d) / (u_m+u_n);
       
        integral = norm_m * norm_n * (outside_wg + inside_wg);
        return
    end
    
end