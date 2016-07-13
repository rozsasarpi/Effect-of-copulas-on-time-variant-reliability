% Sudret et al. (2002). Comparison of methods for computing the probability of failure in time-variant reliability using the outcrossing approach
% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
% modified version of the published example
% simply supported bridge subjected to stochastic load (F) and corrosion


clear all
close all
clc

%==========================
% OPTIONS
%==========================

% time instant for calculation
t = 1;
% limit state to check (1: t; 2: t + delta_t)
als = 1;

%==========================
% ANALYSIS
%==========================

[probdata, analysisopt, gfundata] = ferum_main_timeint;

% This function updates probdata and gfundata before any analysis (must be run only once) 
[probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

% This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable 
probdata.marg = distribution_parameter(probdata.marg); 
% FORM analysis % 

[formresults, probdata] = form(1, probdata, analysisopt, gfundata, [], []);

beta     = formresults.beta


