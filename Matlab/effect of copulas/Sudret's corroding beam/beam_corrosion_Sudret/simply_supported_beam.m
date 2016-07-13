% Sudret et al. (2002). Comparison of methods for computing the probability of failure in time-variant reliability using the outcrossing approach
% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
% simply supported beam subjected to stochastic load (F) and corrosion

% number 2 in variable names indicates the old equation (direct finite difference)

clearvars
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

% time instant for calculation
t = [0, 2, 5, 10, 15, 20];
% limit state to check (1: t; 2: t + delta_t)
als = [1, 2];

%phi2 = @(x1, x2, corr) 1/(2*pi*sqrt(1-corr^2)) * exp(-1/(2*(1-corr^2)) * (x1.^2 - 2*corr*x1.*x2 + x2.^2));
%==========================
% ANALYSIS
%==========================
psi     = @(x) normpdf(x) - x.*normcdf(-x);

N       = length(t);
beta1   = zeros(N, 1);
Pf_t    = zeros(N, 1);
Pf_t2   = zeros(N, 1);
nu_p    = zeros(N, 1);
nu_p2   = zeros(N, 1);
for j = 1:N

    alpha = zeros(5, 2);
    beta  = zeros(2, 1);
    for i = 1:2

        [probdata, analysisopt, gfundata] = ferum_main(t(j), als(i), ro_F);

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
    nu_p2(j)    = mvncdf([beta(1),-beta(2)],[0,0],[1, ro_appr; ro_appr, 1])/delta_t;

    dbeta_dt    = norm(beta(2)-beta(1))/delta_t; 
    dalpha_dt   = norm(alpha(:,2)-alpha(:,1))/delta_t;
    nu_p(j)     = dalpha_dt*normpdf(beta(1))*psi(dbeta_dt/dalpha_dt);
    beta1(j)    = beta(1);
    
    % mean number of outcrossings
    if t(j) == 0
        mean_nu_p = 0;
        mean_nu_p2= 0;
    else
        mean_nu_p   = trapz(t(1:j), nu_p(1:j));
        mean_nu_p2  = trapz(t(1:j), nu_p2(1:j));
    end
    
    Pf_t(j)     = normcdf(-beta1(1)) + mean_nu_p;
    Pf_t2(j)    = normcdf(-beta1(1)) + mean_nu_p2;

end

[nu_p, nu_p2]
beta_t = -norminv(Pf_t);
beta_t2= -norminv(Pf_t2);

%==========================
% PLOT
%==========================

f_diag = figure('Position',[200, 200, 1200, 500]);
subplot(1,2,1)
plot(t, nu_p, 'o-', 'Color', 'green')
hold on
plot(t, nu_p2, 'o-', 'Color', 'red')
legend('new Eq.', 'old Eq.')

title('Outcrossing rate in time')
ylim([0, 0.02])
set(gca, 'YTick', [0:0.002:0.02])
set(gca, 'XTick', [0:5:20])
ylabel('\nu^+(t)')
xlabel('t [year]')
grid on

subplot(1,2,2)
plot(t, beta_t, 'o-', 'Color', 'green')
hold on
plot(t, beta_t2, 'o-', 'Color', 'red')
xa = [0:2:20];
ya = [4.53, 2.75, 2.47, 2.25, 2.07, 1.82, 1.75, 1.6, 1.44, 1.30, 1.12];
plot(xa, ya, 'x-')

title('Reliability index (upper bound) in time')
legend('PHI2 method (new Eq.)','PHI2 method (old Eq.)', 'by eye from ref. paper')
ylim([0.5, 5])
set(gca, 'YTick', [0.5:0.5:5])
set(gca, 'XTick', [0:5:20])
ylabel('\beta_g_e_n(t)')
xlabel('t [year]')
grid on

