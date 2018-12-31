function a = A(DWG, betas, z)
fprintf('z = %f\n', z)

dim = size(betas, 1) - 1;
a = zeros(dim);

for m = 1 : dim
    for n = 1 : m
        a(n,m) = Anm(DWG, betas, z, n, m);
        a(m,n) = -real(a(n,m)) + imag(a(n,m)) * 1j;
    end
end
end