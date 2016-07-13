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

clear all
close all
clc

%--------------------------------------------------------------------------
% OPTIONS
%--------------------------------------------------------------------------

% COPULA
% copulatype = {'Gaussian', 't', 'Clayton', 'Gumbel'};
% copulatype = {'Gaussian', 't', 'Gumbel'};
copulatype = {'Gaussian'};
% copulatype = {'t'};
% copulatype = {'Clayton'};
% copulatype = {'CClayton'};
% copulatype = {'Gumbel'};
% copulatype = {'HR'};
% copulatype = {'Plackett'};
% 
% copulatype = {'Frank'}; % NOT WORKING

% AUTOCORRELATION
Corr.type   = 'gauss';
Corr.pow    = 2;

% survival probability of a component
Ps_comp     = 1 - 1e-6;

% correlation 'length' for S process
tau_Fv      = 10/365;

% increments
% frac        = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
frac        = logspace(-4,-1,20);

%--------------------------------------------------------------------------
% ANALYSIS/CALCULATION
%--------------------------------------------------------------------------
mm  = numel(copulatype);
nn  = numel(frac);

Pf      = nan(mm, nn);
NU_p    = nan(mm, nn);
NU_pp   = NU_p;
PF_SYS  = nan(mm, nn);
% loop over copula types
for jj = 1:mm
    
    copula = copulatype{jj};
    
    % loop over correlation increments
    for ii = 1:nn
        
        tau_F   = tau_Fv;
        Corr.length = tau_F;
        % time increment
        delta_t = tau_F*frac(ii);
        
        % ktau correlation between 'adjactent' F realizations
%         k_tau   = 2/pi*asin(exp(-(delta_t/tau_F).^2)); % to get back the typical correlation function for pearson rho
                CM = cov_matrix([0, delta_t], Corr);
                k_tau = 2/pi*asin(CM(1,2));
        
        % Pearson correlation for Gauss and t copulas
        ro_F    = sin(k_tau*pi/2);
        
        %------------------------------------------------------------------
        % DIRECT INTEGRATION - INITIALIZATION
        %------------------------------------------------------------------
        % copula
        switch copula
            case 'Gaussian'
                theta = ro_F;
            case 't' % dof = 2
                theta = ro_F;
            case 'Clayton'
                theta = 2*k_tau/(1-k_tau);
%                 theta = 1/delta_t;
            case 'CClayton'
                theta = 2*k_tau/(1-k_tau); %???
            case 'Frank'
                theta = fzero(@(alpha) copulastat('Frank', alpha, 'type', 'Kendall') - k_tau, k_tau);
            case 'Plackett'
                theta = fzero(@(alpha) copulastat('Plackett', alpha, 'type', 'Kendall') - k_tau, k_tau);
            case 'Gumbel'
                theta = 1/(1-k_tau);
            case 'HR'
%                 theta = hr_ktau2delta(k_tau);
        end
        
        %------------------------------------------------------------------
        % DIRECT INTEGRATION - CALCULATION - BUILT-IN COPULACDF
        %------------------------------------------------------------------
        
        if strcmp(copula, 't')
            Pf_sys  = Ps_comp - copulacdf(copula, [Ps_comp, Ps_comp], theta, 2);
        else
            theta
            Pf_sys  = Ps_comp - copulacdf(copula, [Ps_comp, Ps_comp], theta);
%             Pf_sys = Ps_comp - bicclay_copulacdf([Ps_comp, Ps_comp], theta);
%             Pf_sys/Pf_sys2
        end
        
        PF_SYS(jj,ii)   = Pf_sys;
        
        nu_p            = Pf_sys/delta_t;
        NU_p(jj,ii)     = nu_p;
        
        NU_pp(jj,ii)    = exp(log(Pf_sys) - log(delta_t));
        %------------------------------------------------------------------
        %------------------------------------------------------------------
%         beta_f          = -norminv(1-Ps_comp);
%         NU_pp(jj,ii)    = mvnpdf([beta_f, -beta_f], [0, 0], [1, theta; theta, 1])...
%             *2*delta_t/tau_F^2*exp(-(delta_t/tau_F)^2);
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
        Pf0 = 1 - Ps_comp;
        Pf(jj,ii) = Pf0 + nu_p*50;
          
    end
end
copulatype

% for Clayton copula, analytic derivative from Wolfram Alpha
delta_t = tau_F*frac;

% nPF_SYS = bsxfun(@rdivide, PF_SYS, PF_SYS(:,1));
% plot(frac, PF_SYS, '-o')
% semilogx(frac, PF_SYS, '-o')
semilogx(frac, Pf, '-o')
hold on
semilogx(frac, Pf0 + NU_pp*50, '--x')
% hold on
% figure
% semilogx(frac, NU_p2, '-o')
% % loglog(frac, PF_SYS, '-o')
% legend(copulatype, 'Location', 'best')
% ylabel('P_{f,sys}(frac)')
% xlabel('frac = \Delta t/corr length')
% grid on

% figure
% nPf = bsxfun(@rdivide, Pf, Pf(:,1));
% semilogx(frac, nPf, '-o')
% legend(copulatype, 'Location', 'best')
% ylabel('P_f(frac)/P_f(min(frac))')
% xlabel('frac = \Delta t/corr length')
% grid on