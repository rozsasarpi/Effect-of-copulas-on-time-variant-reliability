clear all
close all
clc

mu              = 1;                     % mean value
sigma           = 0.2*mu;                % standard deviation
corr_length     = 10;                     % correlation length
pow             = 2;                     % exponent to calculate cov_mx

L_meshed        = 1:0.1:10;
N               = 10;                   % simulation number

[cov_matrix, L_midpoint] = element_cov_matrix(L_meshed, corr_length, pow);

r = mvnrnd(ones(size(cov_matrix,1),1)*mu, cov_matrix*sigma^2, N)';

plot(repmat(L_midpoint,1,N), r)

