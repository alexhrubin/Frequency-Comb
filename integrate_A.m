%z_prime is <= waveguide's length
function integratedA = integrate_A(WG, z_prime)
step_size = 1 / WG.steps;

total_A = zeros(2);

    for i = 1 : z_prime*WG.steps+1
        a = A(WG, fix((i-1)/WG.steps));
        total_A = total_A + a;
    end
    
integratedA = total_A * step_size;
end
