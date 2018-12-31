function [z,c] = solve_c(DWG, z_max, c0)

zq = 0:0.5:z_max;
betas = all_betas(DWG, zq);
betas = cat(1, betas(1:length(c0), :), betas(end, :));

zspan = [0 z_max];

opts = odeset('AbsTol',1e-3);

tic
[z, c] = ode45(@(z,c) dcdz(z,c,DWG,betas), zspan, c0, opts);
fprintf('time to solve = %f\n', toc)

figure
xlabel('z along waveguide (microns)')
ylabel('absolute value of mode coefficient')
hold on

legendCell = {};

for order = 1 : size(c, 2)
    plot(z, abs(c(:, order)))
    legendCell{order} = num2str(order);
end

%plot(z, arrayfun(@(z) WG.half_width(z), z));
%legendCell{end+1} = 'waveguide half-width';

title('abs(c_i)')
legend(legendCell, 'Location', 'northwest')
grid on

end