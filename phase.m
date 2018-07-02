% Calculates the phase portion of the differential equation
% exp( i * integral(beta_m(z') - beta_n(z')) dz')
% We want ~1
function p = phase (WG, m, n, z)
    
    bd = WG.betas(m,:) - WG.betas(n,:);

    step_size = 0.05; %this is set in Waveguide.m as a result of using 20 points/um
    
    integral = sum(bd) * step_size;
    p = exp(1j * integral);

end