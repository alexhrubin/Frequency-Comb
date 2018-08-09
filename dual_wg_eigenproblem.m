function err = dual_wg_eigenproblem(beta, k, eps, a1, a2, sep)

eta = sqrt(beta^2 - k^2);
gam = sqrt(eps*k^2 - beta^2);

left_side = exp(-2*eta*sep) * k^2*(eps-1)^2 * sin(gam*a1)*sin(gam*a2) + (2*beta^2 +k^2*(eps-1))*sin(gam*a1);

right_side = (2*eta*gam*cos(gam*a1) + (eta^2-gam^2)*sin(gam*a1)) * (2*eta*gam*cos(gam*a2) + (eta^2-gam^2)*sin(gam*a2));

err = left_side - right_side;
end

%right_side = (2*eta*gam*cos(gam*a1) + (2*beta^2 + k^2*(eps-1))*sin(gam*a1)) * (2*eta*gam*cos(gam*a2) + (2*beta^2 + k^2*(eps-1))*sin(gam*a2));
