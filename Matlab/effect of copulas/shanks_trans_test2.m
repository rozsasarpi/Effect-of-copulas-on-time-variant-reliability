clear variables
close all
clc
format longE

% fun             = @(x) sin(x)./x;

tau_Fv = [0.1, 1, 10, 100, 1000]/365;
fun = @(x) direct_int_ktau_fun(delta_t, tau_Fv(1));

% Calculate limit with Shanks transformation
[SA, A, step]   = shanks_transform(fun, 2);

idx             = ~isnan(SA(:,end));
tmp             = SA(idx,end);

step(end)
slim            = tmp(end)