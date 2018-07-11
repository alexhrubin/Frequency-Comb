% This is the integral of dS E_m*E_n
function integral = overlap_integral(WG, z, n, m)
    d = WG.width(z);
    beta_n = WG.betas(n, z*WG.steps+1);
    beta_m = WG.betas(m, z*WG.steps+1);
    
    u_n = sqrt( (WG.n*WG.k)^2 - beta_n^2);
    w_n = sqrt( (beta_n^2 - WG.k^2));
    
    u_m = sqrt( (WG.n*WG.k)^2 - beta_m^2);
    w_m = sqrt( (beta_m^2 - WG.k^2));
    
    
    
    if (mod(n,2) == 0) && (mod(m,2) ==1) %symm * asymm
        integral = 2*sin(u_n*d)*cos(u_m*d) / (w_n+w_m);
        return
    
    elseif (mod(m,2) == 0) && (mod(n,2) ==1) %symm * asymm
        integral = 2*cos(u_n*d)*sin(u_m*d) / (w_n+w_m);
        return
     
        
    elseif (mod(n,2)==1) && (mod(m,2)==1) %symmetric modes
        outside_wg = 2*cos(u_n*d)*cos(u_m*d) / (w_n+w_m);
        C_nm = cos(u_n*d)*exp(w_n*d) * cos(u_m*d)*exp(w_m*d);
        inside_wg = C_nm * (2*u_m*sin(u_n*d)*cos(u_m*d) - 2*u_n*cos(u_n*d)*sin(u_m*d)) / (u_n^2 - u_m^2);
        
        integral = outside_wg + inside_wg;
        return
        
 
    elseif (mod(n,2)==0) && (mod(m,2)==0) %anti-symmetric modes
        outside_wg = 2*sin(u_n*d)*sin(u_m*d) / (w_n+w_m);
        C_nm = sin(u_n*d)*exp(w_n*d) * sin(u_m*d)*exp(w_m*d);
        inside_wg = C_nm * (2*u_n*sin(u_n*d)*cos(u_m*d) - 2*u_m*cos(u_n*d)*sin(u_m*d)) / (u_n^2 - u_m^2);
       
        integral = outside_wg + inside_wg;
        return
    end
    
end