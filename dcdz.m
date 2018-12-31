function d = dcdz(z, c, DWG, betas)

a = A(DWG, betas, z);
d = a * c;
end
