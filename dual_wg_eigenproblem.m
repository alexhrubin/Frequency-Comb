function err = dual_wg_eigenproblem(beta, k, eps, a1, a2, sep)
% here, sep is the distance between the inner edges of the waveguides

gam = sqrt(eps*k^2 - beta^2);
eta = sqrt(beta^2 - k^2);

left_side = exp(-2*eta*sep) * (eta^2+gam^2)^2 * sin(gam*a1)*sin(gam*a2);

right_side = (2*eta*gam*cos(gam*a1) + (eta^2-gam^2)*sin(gam*a1)) * (2*eta*gam*cos(gam*a2) + (eta^2-gam^2)*sin(gam*a2));

err = left_side - right_side;
end