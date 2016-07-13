% Check the effect of finite difference step on approximation of
% derivative!
%
% Time-variant reliability analysis with continuous stochastic process
% investigation of the effect of copula assumption
%
%
% g = R - S(t)
%
% R    - constant
% S(t) - continuous stochastic process

clearvars
close all
clc

%--------------------------------------------------------------------------
% OPTIONS
%--------------------------------------------------------------------------

mp.Digits(100);
% COPULA
% copulatype = {'Gaussian', 't', 'Clayton', 'Gumbel'};
% copulatype = {'Gaussian', 't', 'Gumbel'};
% copulatype = {'Gaussian'};
% copulatype = {'t'};
% copulatype = {'Clayton'};
copulatype = {'rotClayton'}; % 180° rotated Clayton copula
% copulatype = {'Gumbel'};
% copulatype = {'rotGumbel'}; % 180° rotated Gumbel copula
% copulatype = {'HR'};
% copulatype = {'Plackett'};
% copulatype = {'Frank'}; 
% copulatype = {'Joe'};

% AUTOCORRELATION
% Corr.type   = 'gauss';
% Corr.pow    = 2;

% survival probability of a component
Ps_comp     = mp('1 - 1e-10');

% correlation 'length' for S process
tau_Fv      = mp('1/365');
% tau_Fv      = mp('[0.1, 1, 10, 100, 1000]')/mp('365');

% increments
% frac        = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
frac        = logspace(mp('-5'),mp('-1'),20);
% frac        = mp('1e-4');

%--------------------------------------------------------------------------
% ANALYSIS/CALCULATION
%--------------------------------------------------------------------------
nn          = numel(frac);

copula      = copulatype{1};

mm          = length(tau_Fv);
for jj = 1:mm
    tau_F       = tau_Fv(jj)
    
    Pf          = nan(1, nn);
    Pf_l        = nan(1, nn);
    Pf_u        = nan(1, nn);
    NU_p        = nan(1, nn);
    NU_p_l      = nan(1, nn);
    NU_p_u      = nan(1, nn);

    NU_p_mp     = mp(zeros(1,nn));
    PF_SYS      = nan(1, nn);
% loop over correlation increments
for ii = 1:nn
    
    Corr.length = tau_F;
    % time increment
    delta_t     = tau_F*frac(ii);
    
    % ktau correlation between 'adjactent' F realizations
    k_tau       = 2/pi*asin(exp(-(delta_t/tau_F).^2)); % to get back the typical correlation function for pearson rho
%     k_tau       = 2/pi*asin((1+(delta_t/tau_F).^2).^(-2));
    
    % Pearson correlation for Gauss and t copulas
    ro_F        = sin(k_tau*pi/2);
    u           = [Ps_comp, Ps_comp];
    %------------------------------------------------------------------
    % DIRECT INTEGRATION - INITIALIZATION
    %------------------------------------------------------------------
    % copula
    switch copula
        case 'Gaussian'
            theta       = ro_F;
            Pf_sys_mp   = Ps_comp - binorm_copulacdf_mp(u, theta);
%             binorm_copulacdf_mp(u, theta)
        case 't' % dof = 2
            theta       = ro_F;
            Pf_sys_mp   = Ps_comp - bit_copulacdf_mp(u, theta, 2);
        case 'Clayton'
            theta       = 2*k_tau/(1-k_tau);
            Pf_sys_mp   = Ps_comp - biclay_copulacdf_mp(u, theta);
            
%             biclay_copulacdf_mp(u, theta)
%             copulacdf(copula, double(u), double(theta))
        case 'rotClayton'
            theta       = 2*k_tau/(1-k_tau); %???
            Pf_sys_mp   = Ps_comp - bicclay_copulacdf_mp(u, theta);
        case 'Frank'
%             theta = fzero(@(alpha) copulastat('Frank', alpha, 'type', 'Kendall') - double(k_tau), double(k_tau));
            % costly if high precision is requested!
            theta       = fzero(@(alpha) 1 + 4 .* (debye_mp(alpha,1)-1) ./ alpha - k_tau, k_tau)
            Pf_sys_mp   = Ps_comp - bifrank_copulacdf_mp(u, theta);
        case 'Plackett'
            % costly if high precision is requested!
            theta = theta_estimationPlackett_Kendall(k_tau)
            Pf_sys_mp   = Ps_comp - biplackett_copulacdf_mp(u, theta);
        case 'Gumbel'
            theta       = 1/(1-k_tau);
            Pf_sys_mp   = Ps_comp - bigumb_copulacdf_mp(u, theta);
        case 'rotGumbel'
            theta       = 1/(1-k_tau);
            Pf_sys_mp   = Ps_comp - (u(1) + u(2) -1 + bigumb_copulacdf_mp([1-u(1), 1-u(2)], theta));    
        case 'HR'
            %                 theta = hr_ktau2delta(k_tau);
        case 'Joe'
            theta       = tau2par_joe_mp(k_tau);
            Pf_sys_mp   = Ps_comp - bijoe_copulacdf_mp(u, theta);
    end
    
    %------------------------------------------------------------------
    % DIRECT INTEGRATION - CALCULATION - BUILT-IN COPULACDF
    %------------------------------------------------------------------
    
    if strcmp(copula, 't')
%         Pf_sys  = Ps_comp - copulacdf(copula, double(u), double(theta), 2);
        Pf_sys = NaN;
    else
        Pf_sys = NaN;
%         Pf_sys  = Ps_comp - copulacdf(copula, double(u), double(theta));
%         Pf_sys_mp = Ps_comp - binorm_copulacdf_mp([Ps_comp, Ps_comp], theta);
        %             Pf_sys = Ps_comp - bicclay_copulacdf([Ps_comp, Ps_comp], theta);
    end
    
    % Copula bounds
    Pf_sys_u    = Ps_comp - bifh_bounds(u, 'lower');
    Pf_sys_l    = Ps_comp - bifh_bounds(u, 'upper');
    
    PF_SYS(ii)  = Pf_sys;
    
    nu_p        = Pf_sys/delta_t;
    NU_p(ii)    = nu_p;
    
    NU_p_mp(ii) = Pf_sys_mp/delta_t
    
    NU_p_l(ii)  = Pf_sys_l/delta_t;
    NU_p_u(ii)  = Pf_sys_u/delta_t;
    
    Pf0         = 1 - Ps_comp;
    Pf(ii)      = Pf0 + nu_p*50;
  
end


% WARNING
Pf(1)
% Pf(12)
Pf_mp = Pf0 + NU_p_mp(1)*50
% for Clayton copula, analytic derivative from Wolfram Alpha
delta_t = tau_F*frac;

figure
semilogx(frac, NU_p, '-o')
hold on
semilogx(frac, NU_p_mp, '--x')
% semilogx(frac, NU_p_l, '.')
% semilogx(frac, NU_p_u, '.')

hl = legend('double', 'multiprecision', 'Location', 'best');
hl.Interpreter = 'LaTeX';
xlabel('$\Delta_\mathrm{t}/\tau_\mathrm{F}$')
ylabel('$\nu^+$')
title(copula, 'Interpreter', 'LaTeX')

set(gca, 'TickLabelInterpreter', 'LaTeX')

end

