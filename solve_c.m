function solve_c(WG, step, zmax, c0)

zq = 0:step:zmax;
betas = all_betas(WG, zq);

zspan = [0 zmax];

[z, c] = ode45(@(z,c) dcdz(z,c,WG,betas,1,3), zspan, c0);

plot(z, abs(c(:,1)));
title c_1

figure
plot(z, abs(c(:,3)));
title c_3
end