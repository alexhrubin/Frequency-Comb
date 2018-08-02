%takes the E, H fields from the maxwellFDFD solver and extracts the
% mode coefficients
function ci = get_coeffs(i, E, H, extra)

dimensions = extra.grid3d.L;
dom_width = dimensions(2);

pixels = extra.grid3d.N;
across_wg = linspace(-dimensions(2)/2, dimensions(2)/2, pixels(2));

WG = Waveguide(sqrt(12), 2*pi/1.55, 0.25, 0.04);

Ex = E{3};
Ex = Ex(round(pixels(1)/2), :);

Hy = H{2};
Hy = Hy(round(pixels(1)/2), :);

E_raw_x = linspace(-dom_width/2, dom_width/2, pixels);
H_raw_x = linspace(-dom_width/2, dom_width/2, pixels);

Ex = interp1(E_raw_x, Ex, across_wg);
Hy = -interp1(H_raw_x, Hy, across_wg);

% hold on
% plot(across_wg, Ex)
% plot(across_wg, Hy)
ci = coeff_integral(WG, dimensions(1), i, Ex, Hy, across_wg);

end