function a = A(DWG, betas, z, starting_guess)
fprintf('z = %f\n', z)

dim = size(betas, 1) - 1;
a = zeros(dim);

for m = 1 : dim
    for n = 1 : m
        a(n,m) = Anm(DWG, betas, z, n, m, starting_guess);
        a(m,n) = -real(a(n,m)) + imag(a(n,m)) * 1j;
    end
end
end

