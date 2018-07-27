function anm = Anm(WG, betas, z, n, m)
disp(z)
if (n == m)
    anm = 0;
    return
end

local_betas = WG.getbeta(z);

if (max(n,m) > length(local_betas))
    anm = 0;
    return
end

numerator = WG.dedz() * WG.k * overlap(WG, z, n, m);
denominator = local_betas(m) - local_betas(n); 

if (z < 0.1)
    z_fine = linspace(0, z, 10);
else
    z_fine = 0:0.1:z;
end

if (z ~= 0)
    zq = betas(end, :);
    beta_diff = interp1(zq, betas(m, :), z_fine) - interp1(zq, betas(n, :), z_fine);
    integral = trapz(z_fine, beta_diff);
    phase = exp(1j * integral);
else
    phase = 1;
end

anm = phase * numerator / denominator;