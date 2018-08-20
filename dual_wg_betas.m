function b= dual_wg_betas(a1, a2, d)

wg1 = Waveguide(12, 2*pi/1.55, @(z) a1/2, @(z) -a1/2);
wg2 = Waveguide(12, 2*pi/1.55, @(z) a2/2, @(z) -a2/2);

naked_wg_betas = sort([wg1.getbeta_position(0), wg2.getbeta_position(0)]);

options = optimset('Display', 'off');
f = @(x) dual_wg_eigenproblem(x, 2*pi/1.55, 12, a1, a2, d);

b = [];
for i = naked_wg_betas
    b = [b, fsolve(f, i, options)];
end
disp(naked_wg_betas - b)
end