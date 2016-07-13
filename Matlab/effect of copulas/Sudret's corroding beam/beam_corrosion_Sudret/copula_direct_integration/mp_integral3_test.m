clearvars
close all
clc

format long
mp.Digits(34);

fun     = @(x,y) y.*sin(x)+x.*cos(y);
ymax    = @(x) x;
tic
integral2(fun, pi,2*pi,0,ymax)
% integral2(fun, mp('pi'),2*pi,0,ymax)
toc

% tic
% dblquad(fun, mp('pi'),2*mp('pi'),0,mp('pi'))
% toc



% % basic example from help
% %% EXAMPLE 1
% fun = @(x,y,z) y.*sin(x)+z.*cos(x)
% 
% % built-in
% tic
% q1 = integral3(fun,0,pi,0,1,-1,1, 'AbsTol',1e-16, 'RelTol',1e-12)
% toc
% 
% % multiprecision
% tic
% q2 = triplequad(fun,0,mp('pi'),0,1,-1,1)
% toc
% 
% %% EXAMPLE 2
% fun = @(x,y,z) x.*cos(y) + x.^2.*cos(z);
% 
% xmin = -1;
% xmax = 1;
% ymin = @(x)     -1*sqrt(1 - x.^2);
% ymax = @(x)      sqrt(1 - x.^2);
% zmin = @(x,y)   -1*sqrt(1 - x.^2 - y.^2);
% zmax = @(x,y)    sqrt(1 - x.^2 - y.^2);
% 
% % built-in
% tic
% q1 = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax,'Method','tiled')
% toc
% 
% % multiprecision
% tic
% q2 = triplequad(fun, -mp('1'),mp('1'),ymin,ymax,zmin,zmax)
% toc