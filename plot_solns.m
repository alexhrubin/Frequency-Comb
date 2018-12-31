function plot_solns(V)
urange = 0:0.01:V;

figure
hold on
grid on
plot(urange, arrayfun(@(u) sqrt(V^2 - u^2), urange))
plot(urange, arrayfun(@(u) u*tan(u), urange))
plot(urange, arrayfun(@(u) -u*cot(u), urange))
ylim([0 V])
hold off
end