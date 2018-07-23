%generate a set of A(z) matrices at the specified query points zq
function As = set_of_A(WG, zq)

betas = all_betas(WG, zq);

As = {};

for i = 2 : length(zq)
    As{i} = A(WG, betas, i);
end
    