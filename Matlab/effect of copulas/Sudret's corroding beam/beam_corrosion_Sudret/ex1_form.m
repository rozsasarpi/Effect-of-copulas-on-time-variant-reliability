% Sudret et al. (2002). Comparison of methods for computing the probability of failure in time-variant reliability using the outcrossing approach
% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
% simply supported bridge subjected to stochastic load (F) and corrosion

% number 2 in variable names indicates the old equation (direct finite difference)

clear all
close all
clc

%==========================
% OPTIONS
%==========================

% correlation 'length' for F process
tau_F = 10/365;
% time increment
delta_t = tau_F*0.1;
% correlation between 'adjactent' F realizations
k_tau = exp(-(delta_t/tau_F));
ro_F = sin(k_tau*pi/2);

% limit state to check (1: t; 2: t + delta_t)
als = [1, 2];

%phi2 = @(x1, x2, corr) 1/(2*pi*sqrt(1-corr^2)) * exp(-1/(2*(1-corr^2)) * (x1.^2 - 2*corr*x1.*x2 + x2.^2));
%==========================
% ANALYSIS
%==========================
psi     = @(x) normpdf(x) - x.*normcdf(-x);


alpha = zeros(2, 2);
beta  = zeros(2, 1);
for i = 1:2
    
    [probdata, analysisopt, gfundata] = ferum_ex1_form(als(i), ro_F);
    
    % This function updates probdata and gfundata before any analysis (must be run only once)
    [probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);
    
    % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
    probdata.marg = distribution_parameter(probdata.marg);
    % FORM analysis %
    
    [formresults, probdata] = form(1, probdata, analysisopt, gfundata, [], []);
    
    beta(i)     = formresults.beta;
    alpha(:,i)  = formresults.alpha;
    
end

ro_appr     = -alpha(:,1)'*alpha(:,2);
nu_p2       = mvncdf([beta(1),-beta(2)],[0,0],[1, ro_appr; ro_appr, 1])/delta_t;

dbeta_dt    = norm(beta(2)-beta(1))/delta_t;
dalpha_dt   = norm(alpha(:,2)-alpha(:,1))/delta_t;
nu_p        = dalpha_dt*normpdf(beta(1))*psi(dbeta_dt/dalpha_dt);
beta1       = beta(1);


[nu_p, nu_p2]




