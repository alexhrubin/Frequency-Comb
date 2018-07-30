function ci = get_coeffs(E, H, solveinfo, )

dimensions = solveinfo.grid3d.L;
dom_width = dimensions(2);

pixels = solveinfo.grid3d.Ncell{2};
across_wg = linspace(-dimensions(2)/2, dimensions(2)/2, pixels);

WG = Waveguide(sqrt(12), 2*pi/1.55, 0.25, 0.02);

E_raw = abs(E{3}.array(end-10, :, 2));
H_raw = abs(H{1}.array(end-10, :, 2));

E_raw_x = linspace(-dom_width/2, dom_width/2, size(E_raw, 2));
H_raw_x = linspace(-dom_width/2, dom_width/2, size(H_raw, 2));

Ex = interp1(E_raw_x, E_raw, across_wg);
Hy = interp1(H_raw_x, H_raw, across_wg);

c_order = coeff_integral(WG, dimensions(1), 1, Ex, Hy, across_wg)

end