function d = dcdz(z, c, WG, betas, max_modes)

a = A(WG, betas, z, max_modes);
d = a * c;
end
