function err=asym_modes(U, V);

%dispersion equations
W=sqrt(V^2-U^2);

%Mode equations
err=abs(-U*cot(U)-W);
