function r = rusting_beam_resistance(fy, b0, h0, r, t)

% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
%[m], [year], [N]

% corrosion rate
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
kappa   = 0.05*10^-3;   % m/year
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% span
L       = 5;
% unit weight of the beam [N/m^3]
ro      = 78.5*10^3;

% reduction of dimensions
b       = b0 - 2*kappa*t;
h       = h0 - 2*kappa*t;

r       = (b.*h.^2.*fy./4 - ro.*b0.*h0.*L.^2/8) - r;

end


