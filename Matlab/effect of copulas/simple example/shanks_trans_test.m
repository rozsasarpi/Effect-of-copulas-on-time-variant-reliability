clear variables
close all
clc
format longE

% fun             = @(x) sin(x)./x;

tau_Fv  = [0.1, 1, 10, 100, 1000]/365;
tau_F   = tau_Fv(3);
fun     = @(dt) direct_int_ktau_fun(dt, tau_F);

fun(0.1*tau_F)
%
% Richardson extrapolation
k       = 10;
ini     = 2;
kk      = ini:k;
seq     = 2.^(-kk).';
A       = arrayfun(fun, seq);

plot(seq/tau_F, A/A(end), 'o-')

% fit polynomial
p       = polyfit(seq, A, k-2);
A
nu_p    = polyval(p, 0)



% Calculate limit with Shanks transformation
% [SA, A, step]   = shanks_transform(fun, 2);
% 
% idx             = ~isnan(SA(:,end));
% tmp             = SA(idx,end);
% 
% step(end)
% slim            = tmp(end)


% nu_p*50