% Calculates the phase portion of the differential equation
% exp( i * integral(beta_m(z') - beta_n(z')) dz')
% We want ~1
function p = phase (WG, z, n, m)
    
    bd = WG.betas(m, 1:z*WG.steps+1) - WG.betas(n, 1:z*WG.steps+1);
    
    step_size = 1 / WG.steps; %this is set in Waveguide.m as a result of using 5 points/um
    
    integral = sum(bd) * step_size;
    p = exp(1j * integral);

end