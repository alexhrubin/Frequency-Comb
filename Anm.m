function anm = Anm(WG, betas, z, n, m)
if (z == 0)
    anm = 0;
    return
end

local_betas = WG.getbeta(z);

numerator = WG.dedz() * WG.k * overlap_numeric(WG, z, n, m);
denominator = local_betas(m) - local_betas(n);

z_to_use = betas(end, :) <= z;
zq = betas(end, z_to_use);

beta_diff = betas(m, z_to_use) - betas(n, z_to_use);
integral = trapz(zq, beta_diff);
phase = exp(1j * integral);

anm = phase * numerator / denominator;