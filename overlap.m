function overlap = overlap(WG, z, n, m)

mode_m = WG.eigenmode_function(z, m);
mode_n = WG.eigenmode_function(z, n);

top = WG.upper_edge_func(z);
bottom = WG.lower_edge_func(z);

upper_slope = WG.upper_edge_dz(z);
lower_slope = WG.lower_edge_dz(z);

mntop = mode_n(top);
mmtop = mode_m(top);

overlap = (WG.eps - 1) * (upper_slope * mode_m(top)*mode_n(top) - lower_slope * mode_m(bottom)*mode_n(bottom));
end
