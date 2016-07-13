clear all
close all
clc


%% TWO-VARIABLE PROBLEM
% g = R - E

fR      = @(r) normpdf(r,4,1);
fE      = @(e) normpdf(e,0,1);
fRE     = @(r,e) fR(r).*fE(e);
emin    = @(r) r;

% exact solution
normcdf(-4/sqrt(2))

tic
% solution with integral2
integral2(fRE, -2, 10, emin, 6)
toc

tic
% solution with nested integral, 2-layer
integral(@(r) innerE(r, fRE), -2, 10, 'ArrayValued', true)
toc

%% THREE-VARIABLE PROBLEM
% g = R - E - S

fS      = @(s) normpdf(s,0,1);
fRES    = @(r,e,s) fR(r).*fE(e).*fS(s);
smin    = @(r,e) r-e;

% exact solution
normcdf(-4/sqrt(3))

tic
% solution with integral3
integral3(fRES, -2, 10, -6, 6, smin, 6)
toc


% way too slow solution > gauss quadrature?
% tic
% % solution with nested integral - 3-layer
% integral(@(r) innerES(r, fRES), -2, 10, 'ArrayValued', true)
% toc
