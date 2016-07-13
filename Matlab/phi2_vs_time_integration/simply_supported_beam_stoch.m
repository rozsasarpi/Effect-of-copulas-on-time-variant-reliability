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

% correlation 'length' for F process
tau_F = 100/365;
% time increment
delta_t = tau_F*0.1;
% correlation between 'adjactent' F realizations
ro_F = exp(-(delta_t/tau_F)^2);

% time instant for calculation
t = [0, 2, 5, 10, 15, 20];
%t = 0:0.2:1;
% limit state to check (1: t; 2: t + delta_t)
als = [1, 2];

%==========================
% ANALYSIS
%==========================
psi     = @(x) normpdf(x) - x.*normcdf(-x);

N       = length(t);
beta1   = zeros(N, 1);
Pf_t    = zeros(N, 1);
nu_p    = zeros(N, 1);
for j = 1:N

    alpha = zeros(5, 2);
    beta  = zeros(2, 1);
    for i = 1:2

        [probdata, analysisopt, gfundata] = ferum_main_stoch(t(j), als(i), ro_F);

        % This function updates probdata and gfundata before any analysis (must be run only once) 
        [probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

        % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable 
        probdata.marg = distribution_parameter(probdata.marg); 
        % FORM analysis % 

        [formresults, probdata] = form(1, probdata, analysisopt, gfundata, [], []);

        beta(i)     = formresults.beta; 
        alpha(:,i)  = formresults.alpha;

    end

    dbeta_dt    = norm(beta(2)-beta(1))/delta_t; 
    dalpha_dt   = norm(alpha(:,2)-alpha(:,1))/delta_t;
    nu_p(j)     = dalpha_dt*normpdf(beta(1))*psi(dbeta_dt/dalpha_dt);
    beta1(j)    = beta(1);
    
    % mean number of outcrossings
    if t(j) == 0
        mean_nu_p = 0;
    else
        mean_nu_p   = trapz(t(1:j), nu_p(1:j));
    end
    
    Pf_t(j)     = normcdf(-beta1(1)) + mean_nu_p;

end


beta_t = -norminv(Pf_t);

%==========================
% PLOT
%==========================

f_diag = figure('Position',[200, 200, 1200, 500]);
subplot(1,2,1)
plot(t, nu_p, 'o-', 'Color', 'red')

title('Outcrossing rate in time')
ylabel('\nu^+(t)')
xlabel('t [year]')
grid on

subplot(1,2,2)
plot(t, beta_t, 'o-', 'Color', 'red')

title('Reliability index (upper bound) in time')

ylabel('\beta_g_e_n(t)')
xlabel('t [year]')
grid on

beta_t(end)
