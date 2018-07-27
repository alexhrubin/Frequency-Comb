function overlap = overlap(WG, z, n, m)

mode_m = WG.eigenmode_function(z, m);
mode_n = WG.eigenmode_function(z, n);
d = WG.width(z);

overlap = (sqrt(WG.n) - 1) * WG.dedz() * 2 * mode_m(d)*mode_n(d);
end
