% This is the integral of dS E_m*E_n
function integral = overlap_integral(profiles, n, m)
    step_size = (profiles(end,end) - profiles(end,1)) / length(profiles);
    integral = sum((profiles(n,:) .* profiles(m,:)) * step_size);
end