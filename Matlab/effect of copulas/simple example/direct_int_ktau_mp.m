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
copulatype = {'Gaussian', 't', 'Gumbel', 'rotClayton', 'rotGumbel'};
% copulatype = {'Gaussian', 't', 'Gumbel'};
% copulatype = {'Gaussian'};
% copulatype = {'t'};
% copulatype = {'Clayton'};
% copulatype = {'rotClayton'}; % 180° rotated Clayton copula
% copulatype = {'Gumbel'};
% copulatype = {'rotGumbel'}; % 180° rotated Gumbel copula
% copulatype = {'HR'};
% copulatype = {'Plackett'};
% copulatype = {'Frank'}; 
% copulatype = {'Joe'};


% Corr.type   = 'gauss';
% Corr.pow    = 2;

% total duration of time [years]
t_tot = 50;

% survival probability of a component
Ps_comp     = mp('1 - 1e-9');

% correlation 'length' for S process
tau_Fv      = mp('1/365');
% tau_Fv      = mp('[0.1, 1, 10, 100, 1000]')/mp('365');

% increment
frac        = mp('1e-5');

%--------------------------------------------------------------------------
% ANALYSIS/CALCULATION
%--------------------------------------------------------------------------

nn          = length(copulatype);
Pf          = nan(1, nn);
NU_p        = nan(1, nn);
% loop over correlation lengths
for ii = 1:nn
    
    copula      = copulatype{ii};
    tau_F       = tau_Fv;%(ii);
    
    Corr.length = tau_F;
    % time increment
    delta_t     = tau_F*frac;
    
    %......................................................................
    % AUTOCORRELATION
    % Gaussian
    ro_F        = exp(-(delta_t/tau_F)^2);
    % Cauchy
%     ro_F      = (1+(delta_t/tau_F)^2).^(-2);
    %......................................................................
    % Kendall tau
    k_tau       = 2*asin(ro_F)/mp('pi');
    
    u           = [Ps_comp, Ps_comp];
    %------------------------------------------------------------------
    % CALCULATION
    %------------------------------------------------------------------
    % copula
    switch copula
        case 'Gaussian'
            theta       = ro_F;
            Pf_sys   = Ps_comp - binorm_copulacdf_mp(u, theta);
        case 't' % dof = 2
            theta       = ro_F;
            Pf_sys   = Ps_comp - bit_copulacdf_mp(u, theta, 2);
        case 'Clayton'
            theta       = 2*k_tau/(1-k_tau);
            Pf_sys   = Ps_comp - biclay_copulacdf_mp(u, theta);   
        case 'rotClayton'
            theta       = 2*k_tau/(1-k_tau); %???
            Pf_sys   = Ps_comp - bicclay_copulacdf_mp(u, theta);
        case 'Frank'
            % costly if high precision is requested!
            theta       = fzero(@(alpha) 1 + 4 .* (debye_mp(alpha,1)-1) ./ alpha - k_tau, k_tau);
            Pf_sys   = Ps_comp - bifrank_copulacdf_mp(u, theta);
        case 'Plackett'
            % costly if high precision is requested!
            theta = theta_estimationPlackett_Kendall(k_tau);
            Pf_sys   = Ps_comp - biplackett_copulacdf_mp(u, theta);
        case 'Gumbel'
            theta       = 1/(1-k_tau);
            Pf_sys   = Ps_comp - bigumb_copulacdf_mp(u, theta);
        case 'rotGumbel'
            theta       = 1/(1-k_tau);
            Pf_sys   = Ps_comp - (u(1) + u(2) -1 + bigumb_copulacdf_mp([1-u(1), 1-u(2)], theta));    
        case 'HR'
            %                 theta = hr_ktau2delta(k_tau);
        case 'Joe'
            theta    = tau2par_joe_mp(k_tau);
            Pf_sys   = Ps_comp - bijoe_copulacdf_mp(u, theta);
        otherwise
            
    end
    
    nu_p        = Pf_sys/delta_t;
    NU_p(ii)    = nu_p;
     
    Pf0         = 1 - Ps_comp;
    Pf(ii)      = Pf0 + nu_p*50;
end

[NU_p(1:2), 0, NU_p(3:end)]'
[Pf(1:2), 0, Pf(3:end)]'

% xx      = tau_Fv*365;
% p       = polyfit(log10(xx), log10(NU_p),1);
% yy      = polyval(p, log10(xx));
% loglog(tau_Fv*365, NU_p, '-o')
% hold on
% loglog(xx, 10.^yy, '--x')
% xlabel('$\tau_\mathrm{F}/365$')
% ylabel('$\nu^+$')
% set(gca, 'TickLabelInterpreter', 'LaTeX')
% 
% figure
% p       = polyfit(log10(xx), log10(Pf),1);
% yy      = polyval(p, log10(xx));
% loglog(tau_Fv*365, Pf, '-o')
% hold on
% loglog(xx, 10.^yy, '--x')
% xlabel('$\tau_\mathrm{F}/365$')
% ylabel('$P_f$')
% set(gca, 'TickLabelInterpreter', 'LaTeX')


