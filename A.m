function a = A(WG, betas, z, max_modes)
fprintf('z = %f\n', z)

%dim = size(betas, 1) - 1;
dim = max_modes;
a = zeros(dim);

for m = 1 : dim
    for n = 1 : m
        a(n,m) = Anm(WG, betas, z, n, m, max_modes);
        a(m,n) = -real(a(n,m)) + imag(a(n,m)) * 1j;
    end
end
end

