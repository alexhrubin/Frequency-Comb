function bd = betadiff(WG, m, n, z)
    betas = getbeta(WG, z);
    bd = betas(m) - betas(n);
end
