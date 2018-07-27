function d = dcdz(z, c, WG, betas)

a = A(WG, betas, z);
d = a * c;
end
