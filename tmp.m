function y = tmp(beta, a1, a2)

    eps = 12;
    k = 4.053667940115862;
    
    gam = sqrt(eps*k^2 - beta^2);
    eta = sqrt(beta^2 - k^2);
    
    first = (-gam * cos(gam*a1) + eta * sin(gam*a1)) * (eta * cos(gam*a2) - gam * sin(gam*a2));
    second = (-eta * cos(gam*a1) + gam * sin(gam*a1)) * (eta * sin(gam*a2) + gam * cos(gam*a2));
    
    y = first + second;
end