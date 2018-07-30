function overlap = overlap(WG, z, n, m)


mode_m = WG.eigenmode_function(z, m);
mode_n = WG.eigenmode_function(z, n);
d = WG.width(z);


overlap = (WG.n^2 - 1) * WG.dedz(z) * (mode_m(d)*mode_n(d) + mode_m(-d)*mode_n(-d));
end
