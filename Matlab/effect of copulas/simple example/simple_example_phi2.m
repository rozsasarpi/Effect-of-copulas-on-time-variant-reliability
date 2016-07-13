% description..

clear all
close all
clc

%==========================
% OPTIONS
%==========================

% correlation 'length' for F process
tau_F = 1/365;
% time increment
delta_t = tau_F*0.1;
% correlation between 'adjactent' F realizations
ro_F = exp(-(delta_t/tau_F)^2);

% 'resistance'; constant
R = 290.135;%307.97;%324.48;%339.91;%354.45;
% limit state to check (1: t; 2: t + delta_t)
als = [1, 2];

% design life [years]
td = 50;

%phi2 = @(x1, x2, corr) 1/(2*pi*sqrt(1-corr^2)) * exp(-1/(2*(1-corr^2)) * (x1.^2 - 2*corr*x1.*x2 + x2.^2));
%==========================
% ANALYSIS
%==========================
psi     = @(x) normpdf(x) - x.*normcdf(-x);

alpha = zeros(2, 2);
beta  = zeros(2, 1);
for i = 1:2
    
    [probdata, analysisopt, gfundata] = ferum_veri_(R, ro_F, als(i));
    
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

Pf_t        = normcdf(-beta1(1)) + nu_p*td;
Pf_t2       = normcdf(-beta1(1)) + nu_p2*td;


beta_t = -norminv(Pf_t);
beta_t2= -norminv(Pf_t2);

disp('improved equation;    numerical derivation')
[Pf_t, Pf_t2]
[beta_t, beta_t2]
