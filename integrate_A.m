function integratedA = integrate_A(WG, z_prime)
step_size = 0.05;

total_A = zeros(2);

for i = 1 : z_prime * 20
    total_A = total_A + A(WG, i);
    
integratedA = total_A * step_size;
end
