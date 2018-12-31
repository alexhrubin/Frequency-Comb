%find the propagation constants at each query point in zq
function betas = all_betas(DWG, zq, starting_guess)

initial_beta_guess = DWG.naked_wg_betas(zq(1));

beta_cells = cell([1, length(zq)]);
beta_cells{1} = DWG.getbeta_position(zq(1), initial_beta_guess);

for i = 2 : length(zq)
    beta_cells{i} = DWG.getbeta_position(zq(i), beta_cells{i-1});
end


max_betas = max(cellfun(@(item) length(item), beta_cells));


%betas is a zero-padded matrix where each column corresponds to the set of
%propagation constants which exist at the zq_interp points along the
%waveguide
betas = zeros(length(zq), max_betas);
for i = 1:length(zq)
    betas(i, 1:length(beta_cells{i})) = beta_cells{i};
end
betas = transp(betas);
betas = cat(1, betas, zq);


    
