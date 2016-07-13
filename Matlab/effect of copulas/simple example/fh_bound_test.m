clearvars
close all
clc

nn = 10;
u1 = linspace(0,1,nn);
u2 = linspace(0,1,nn);
[U1, U2] = meshgrid(u1, u1);

Cl = reshape(bifh_bounds([U1(:), U2(:)], 'lower'), size(U1));
Cu = reshape(bifh_bounds([U1(:), U2(:)], 'upper'), size(U1));

subplot(1,2,1)
surfc(U1, U2, Cl)
alpha(0.5)
title('lower')

subplot(1,2,2)
surfc(U1, U2, Cu)
colormap('gray')
alpha(0.5)
title('upper')