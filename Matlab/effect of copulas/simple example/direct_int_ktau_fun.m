% Time-variant reliability analysis with continuous stochastic process
% investigation of the effect of copula assumption
%
% g = R - S(t)
%
% R    - constant
% S(t) - continuous stochastic process

function [nu_p, Pf0] = direct_int_ktau_fun(delta_t, tau_Fv)
%--------------------------------------------------------------------------
% OPTIONS
%--------------------------------------------------------------------------

% COPULA
% copulatype = {'Gaussian', 't', 'Clayton', 'Gumbel'};
% copulatype = 'Gaussian';
% copulatype = 't';
copulatype = 'Clayton';
% copulatype = 'Gumbel';

% copulatype = 'Frank'; % NOT WORKING

% AUTOCORRELATION
Corr.type = 'g';
Corr.pow  = 2;

% RESISTANCE
% resistance, to reach the desired Pf0 value, has no effect, just for
% convinience
R = 339.91;%290.135;%640.7;

%--------------------------------------------------------------------------
% ANALYSIS/CALCULATION
%--------------------------------------------------------------------------

copula = copulatype;

tau_F       = tau_Fv;
Corr.length = tau_F;


% ktau correlation between 'adjactent' F realizations
%         k_tau   = 2/pi*asin(exp(-(delta_t/tau_F).^2)); % to get back the typical correlation function for pearson rho
CM      = cov_matrix([0, delta_t], Corr);
k_tau   = 2/pi*asin(CM(1,2));

% Pearson correlation for Gauss and t copulas
ro_F    = sin(k_tau*pi/2);

% marginal properties, have no effect, just for convinience
meanS   = 100;
covS    = 0.4;

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
    case 'Frank'
        theta = fzero(@(alpha) copulastat('Frank', alpha, 'type', 'Kendall') - k_tau, ro_F);
    case 'Plackett'
        % not yet tested
        %fx1x2 = @(x1, x2) theta*(1+(theta+1)*(Fx1(x1)+Fx2(x2)-2*Fx1(x1)*Fx2(x2)))/((1+(theta-1)*(Fx1(x1)+Fx2(x2)))^2-4*Fx1(x1)*Fx2(x2)*theta*(theta-1))^(3/2);
    case 'Gumbel'
        theta = 1/(1-k_tau);
end

%------------------------------------------------------------------
% DIRECT INTEGRATION - CALCULATION - BUILT-IN COPULACDF
%------------------------------------------------------------------
% survival probability of a component
Ps_comp     = normcdf(R, meanS, covS*meanS);
if strcmp(copula, 't')
    Pf_sys  = Ps_comp - copulacdf(copula, [Ps_comp, Ps_comp], theta, 2);
else
    Pf_sys  = Ps_comp - copulacdf(copula, [Ps_comp, Ps_comp], theta);
end

nu_p = Pf_sys/delta_t;

Pf0 = 1 - Ps_comp;
end