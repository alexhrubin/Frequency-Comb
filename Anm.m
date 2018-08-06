function anm = Anm(WG, betas, z, n, m, max_modes)

if (n == m)
    anm = 0;
    return
end

local_betas = WG.getbeta_position(z, max_modes);

if (max(n,m) > length(local_betas))
    anm = 0;
    return
end

if (max(local_betas(m), local_betas(n)) == 0);
    anm = 0;
    return
end

numerator = WG.k * overlap(WG, z, n, m);
denominator = local_betas(m) - local_betas(n);

if (z == 0)
    phase = 1;
else   
    relevant_betas = betas(:, betas(end,:) <= z);
    zero_padding = zeros(1, size(relevant_betas, 1) - length(local_betas) - 1);
    last_col = cat(2, local_betas, zero_padding, z);
    relevant_betas = cat(2, relevant_betas, transp(last_col));
    beta_diff = relevant_betas(m, :) - relevant_betas(n, :);
    integral = trapz(relevant_betas(end, :), beta_diff);

    phase = exp(1j * integral);
end

anm = phase * numerator / denominator;