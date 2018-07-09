% This is the integral of dS E_m*E_n
function integral = overlap_integral(WG, d, n, m, beta_n, beta_m)
    u_n = sqrt( (WG.n*WG.k)^2 - beta_n^2);
    w_n = sqrt( (beta_n^2 - WG.k^2));
    
    u_m = sqrt( (WG.n*WG.k)^2 - beta_n^2);
    w_m = sqrt( (beta_n^2 - WG.k^2));
    
    C_nm = (cos(u_n*d) * exp(w_n*d)) * (cos(u_m*d) * exp(w_m*d));
    
    
    if (mod(n,2) == 0) && (mod(m,2) ==1)
        integral = 0;
        return
        
    elseif (mod(n,2)==1) && (mod(m,2)==1)
        outside_wg = 2*C_nm / (w_n+w_m);
        inside_wg = C_nm * (2*u_m*sin(u_n*d)*cos(u_m*d) - 2*u_n*cos(u_n*d)*sin(u_m*d)) / (u_n^2 - u_m^2);
        
        integral = outside_wg + inside_wg;
        return
        
 
    elseif (mod(n,2)==0) && (mod(m,2)==0)
        outside_wg = 2*C_nm / (w_n+w_m);
        inside_wg = C_nm * (2*u_n*sin(u_n*d)*cos(u_m*d) - 2*u_m*cos(u_n*d)*sin(u_m*d)) / (u_n^2 - u_m^2);
       
        integral = outside_wg + inside_wg;
        return
    end
    
    
   % step_size = (profiles(end,end) - profiles(end,1)) / length(profiles);
   % integral = sum((profiles(n,:) .* profiles(m,:)) * step_size);
end