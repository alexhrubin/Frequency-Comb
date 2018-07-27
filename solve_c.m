function solve_c(WG, z_max, c0)

zq = 0:z_max;
betas = all_betas(WG, zq);

zspan = [0 z_max];

[z, c] = ode45(@(z,c) dcdz(z,c,WG,betas), zspan, c0);

figure
xlabel('z along waveguide (microns)')
ylabel('absolute value of mode coefficient')
hold on

legendCell = {};

for order = 1 : size(c, 2)
    plot(z, abs(c(:, order)))
    legendCell{order} = num2str(order);
end

title('abs(c_i)')
legend(legendCell)


end