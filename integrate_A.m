function integratedA = integrate_A(WG, z_prime)
Z = linspace(0, z_prime, 5);
step_size = z_prime / 5;
   
A_values = arrayfun(@(z) A(WG, z), Z, 'UniformOutput', false);

total_A = zeros(size(A_values{1,1}));
for i = 1:5
    total_A = total_A + A_values{1, i};
    
integratedA = total_A * step_size;
end
