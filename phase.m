% Calculates the phase portion of the differential equation
% exp( i * integral(beta_m(z') - beta_n(z')) dz')
function p = phase (WG, m, n, z)
    num_of_points = 5000;
    Z = linspace(0, z, num_of_points);

    bd = arrayfun(@(zprime) betadiff(WG, m, n, zprime), Z);
    
    integral = sum(bd) * (z/num_of_points);
    p = exp(1j * integral);
end
        
        