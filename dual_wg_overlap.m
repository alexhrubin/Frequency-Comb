function overlap = dual_wg_overlap(DWG, z, n, m)

local_betas = DWG.getbeta_position(z, DWG.naked_wg_betas(z));
beta_m = local_betas(m);
beta_n = local_betas(n);
mode_m = DWG.eigenmode_function(z, beta_m);
mode_n = DWG.eigenmode_function(z, beta_n);

wg1_lower_slope = DWG.wg1.lower_edge_dz(z);
wg1_upper_slope = DWG.wg1.upper_edge_dz(z);
wg2_lower_slope = DWG.wg2.lower_edge_dz(z);
wg2_upper_slope = DWG.wg2.upper_edge_dz(z);

edge1 = DWG.wg1_lower_edge(z);
edge2 = DWG.wg1_upper_edge(z);
edge3 = DWG.wg2_lower_edge(z);
edge4 = DWG.wg2_upper_edge(z);

overlap = (DWG.eps-1) * (wg1_lower_slope*mode_m(edge1)*mode_n(edge1) - wg1_upper_slope*mode_m(edge2)*mode_n(edge2) + wg2_lower_slope*mode_m(edge3)*mode_n(edge3) - wg2_upper_slope*mode_m(edge4)*mode_n(edge4));

%overlap = (DWG.eps - 1) * (upper_slope * mode_m(top)*mode_n(top) - lower_slope * mode_m(bottom)*mode_n(bottom));
end