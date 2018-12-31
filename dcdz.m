function d = dcdz(z, c, DWG, betas, starting_guess)

a = A(DWG, betas, z, starting_guess);
d = a * c;
end
