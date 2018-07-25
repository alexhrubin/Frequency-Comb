wg = Waveguide(sqrt(12), 2*pi/1.55, 0.25, 0.003)

z = 0:0.05:100;

%Let's look at the overlap of modes 1 and 3
%compare the analytic solution with good numerical approximation
overlaps_numerical = arrayfun(@(z) overlap_numeric(wg, z, 1, 3, 0.01), z);
overlaps_analytical = arrayfun(@(z) overlap_integral(wg, z, 1, 3), z);

%they agree, but it just looks like noise...
plot(z, overlaps_numerical)
hold on
plot(z, overlaps_analytical)


%now look at a lower-resolution numerical approximation
overlaps_numerical_low_res = arrayfun(@(z) overlap_numeric(wg, z, 1, 3, 0.05), z)
figure
plot(z, overlaps_numerical_low_res)
% this makes more sense, and we can compare with beta_1 - beta_3 and it is
% convincing: the magnitude of the oscillations decreases as the betas get
% further apart
figure
betas = all_betas(wg, z);
beta_difference = betas(1, :) - betas(3, :);
plot(z, beta_difference)