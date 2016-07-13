% Richardson extrapolation of limit

clear variables
close all
clc

fun     = @(x) sin(x)./x;

lim_x   = 0;
k       = 6;

kk      = 0:k;

% 
seq     = 2.^(-kk);
% seq     = kk.^(-2);

f       = fun(seq);

% fit polynomial
p       = polyfit(seq, f, k-2);

xx      = linspace(lim_x, seq(1),1e3);
pf      = polyval(p, xx);


plim    = polyval(p, lim_x);
1-f(end)
1-plim

min(seq)

plot(seq, f, 'o')
hold on
plot(xx, pf)

