function overlap = overlap(WG, z, n, m)

mode_m = WG.eigenmode_function(z, m);
mode_n = WG.eigenmode_function(z, n);

top = WG.upper_edge(z);
bottom = WG.lower_edge(z);

overlap = (WG.eps - 1) * (WG.dedz_upper(z) * mode_m(top)*mode_n(top) - WG.dedz_lower(z) * mode_m(bottom)*mode_n(bottom));
end
