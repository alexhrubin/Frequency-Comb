function d = dcdz(z, c, WG, betas, n, m)
a13 = Anm(WG, betas, z, n, m);
a31 = Anm(WG, betas, z, m, n);
%-real(a13) + (imag(a13)) * 1j;

d = [a13*c(3);
     0;
     a31*c(1)];
end
