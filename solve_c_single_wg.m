function [z,c] = solve_c(WG, z_max, max_modes, c0)

zq = 0:0.25:z_max;
betas = all_betas(WG, zq, max_modes);

zspan = [0 z_max];

opts = odeset('AbsTol',1e-3);

tic
[z, c] = ode45(@(z,c) dcdz(z,c,WG,betas,max_modes), zspan, c0, opts);
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