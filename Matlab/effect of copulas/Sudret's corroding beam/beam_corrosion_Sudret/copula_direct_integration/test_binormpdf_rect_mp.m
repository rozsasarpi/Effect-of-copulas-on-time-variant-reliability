clearvars
close all
% clc

MU      = mp('[3,2]');
SIGMA   = mp('[1, 0.5; 0.5, 2]');

XL1 = mp('1');%mp('[1, 2; 1, 2]');
XL2 = mp('[2, 2; 2, 2]');
XU1 = mp('[2, 3; 2, 3]');
XU2 = mp('[3, 4; 3, 5]');

rho         = SIGMA(2)/sqrt(SIGMA(1)*SIGMA(4));

tic
binorm_rect_mp(XL1, XL2, XU1, XU2, MU, SIGMA)
toc

tic
bit_rect_mp(XL1, XL2, XU1, XU2, MU, sqrt([SIGMA(1), SIGMA(4)]), rho, 1000)
toc

%--------------------------------------------------------------------------
% Fx1x2   = @(xl, xu)...
%     mvncdf(xu, MU, SIGMA) -...
%     mvncdf([xu(1),xl(2)], MU, SIGMA) -...
%     mvncdf([xl(1),xu(2)], MU, SIGMA) +...
%     mvncdf(xl, MU, SIGMA);
% 
% k = 3;
% Fx1x2([XL1(k), XL2(k)], [XU1(k), XU2(k)])
