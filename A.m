function a = A(WG, betas, z)

dim = size(betas, 1) - 1;
a = zeros(dim);

for m = 1 : dim
    for n = 1 : m
        a(n,m) = Anm(WG, betas, z, n, m);
        a(m,n) = -real(a(n,m)) + imag(a(n,m)) * 1j;
    end
end
