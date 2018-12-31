%find the propagation constants at each query point in zq
function betas = all_betas(WG, zq, max_num_of_modes)

if (nargin == 3)
    beta_cells = arrayfun(@(z) WG.getbeta_position(z, max_num_of_modes), zq, 'UniformOutput', false);
else
    beta_cells = arrayfun(@(z) WG.getbeta_position(z), zq, 'UniformOutput', false);
    max_num_of_modes = max(cellfun(@(item) length(item), beta_cells));
    disp(max_num_of_modes)
end

%betas is a zero-padded matrix where each column corresponds to the set of
%propagation constants which exist at the zq_interp points along the
%waveguide
betas = zeros(length(zq), max_num_of_modes);
for i = 1:length(zq)
    betas(i, 1:length(beta_cells{i})) = beta_cells{i};
end
betas = transp(betas);
betas = cat(1, betas, zq);  


    
