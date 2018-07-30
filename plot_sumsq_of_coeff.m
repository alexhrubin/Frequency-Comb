function plot_sumsq_of_coeff(z, c)

sumsq = [];
for i = 1 : size(c, 1)
    sumsq(i) = sum(abs(c(i,:)) .^ 2);
end

plot(z, sumsq)
end
