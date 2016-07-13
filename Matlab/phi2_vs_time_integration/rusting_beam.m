function g = rusting_beam(F1, F2, fy, b0, h0, t, als)

% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
% simply supported bridge subjected to stochastic load (F) and corrosion
%[m], [year], [N]

% corrosion rate
kappa   = 0;%0.05*10^-3;   % m/year
% span
L       = 5;
% unit weight of the beam [N/m^3]
ro      = 78.5*10^3;

% reduction of dimensions
b       = b0 - 2*kappa*t;
h       = h0 - 2*kappa*t;


switch als
    case 1
        F = F1;
    case 2
        F = F2;
end

% plastic hinge at the middle
g       = b.*h.^2.*fy./4 - ((F.*L)/4 + ro.*b.*h.*L.^2/8);

end


