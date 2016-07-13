clear all
close all
clc

% correlation 'length' for S process
tau_F       = 1/365;
frac        = logspace(-8,-1,20);
delta_t     = tau_F*frac;
u           = 1-1e-8; %Ps_comp

nn          = length(frac);
PSS1        = nan(nn,1);
nu_p2       = nan(nn,1);
for ii = 1:nn
    delta_ti = delta_t(ii);
%     k_tau   = 2/pi*asin(exp(-(delta_ti/tau_F).^2));
    k_tau   = exp(-(delta_ti/tau_F).^2); %!!!
    theta   = 2*k_tau/(1-k_tau);
    
    PSS1(ii) = biclay_copulacdf([u,u], theta);
    
    x = delta_ti;
%     F2(ii) = (2*u^(-theta) - 1)^(-1/theta);
end

Pf_sys1 = u - PSS1;

nu_p1   = Pf_sys1./delta_t';

Pf0 = 1 - u;

Pf1 = Pf0 + nu_p1*50;
% Pf2 = Pf0 + nu_p2*50;

% plot(delta_t, nu_p1)
semilogx(delta_t, nu_p1)
% hold on
% plot(delta_t, nu_p1)
% semilogx(delta_t, nu_p2)
% legend('F1', 'F2')

% Pf1./Pf2

% (2*u^(-theta) - 1)^(-1/theta)
% theta = 2*k_tau/(1-k_tau)
% k_tau = 2/pi*asin(exp(-(x/tau_F)^2))
% exp(-x^2)

% (2 *u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi (1-(2 *asin(exp(-x^2/tau_F^2)))/pi)))-1)^(-(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))/(4 *asin(exp(-x^2/tau_F^2)))) *((-(pi *x *exp(-x^2/tau_F^2) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))/(2 *tau_F^2 *sqrt(1-exp(-(2 *x^2)/tau_F^2)) *asin(exp(-x^2/tau_F^2))^2)-(x *exp(-x^2/tau_F^2))/(tau_F^2 *sqrt(1-exp(-(2 *x^2)/tau_F^2)) *asin(exp(-x^2/tau_F^2)))) *log(2 *u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)))-1)-(pi *log(u) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi) *((8 *x *exp(-x^2/tau_F^2))/(pi *tau_F^2 *sqrt(1-exp(-(2 x^2)/tau_F^2)) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))+(16 *x *exp(-x^2/tau_F^2) *asin(exp(-x^2/tau_F^2)))/(pi^2 *tau_F^2 *sqrt(1-exp(-(2 x^2)/tau_F^2)) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)^2)) *u^(-(4 asin(exp(-x^2/tau_F^2)))/(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))))/(2 *asin(exp(-x^2/tau_F^2)) *(2 *u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)))-1)));


F = @(x) (2*u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)))-1)^(-(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))/(4 *asin(exp(-x^2/tau_F^2)))) *((-(pi *x *exp(-x^2/tau_F^2) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))/(2 *tau_F^2 *sqrt(1-exp(-(2 *x^2)/tau_F^2)) *asin(exp(-x^2/tau_F^2))^2)-(x *exp(-x^2/tau_F^2))/(tau_F^2 *sqrt(1-exp(-(2 *x^2)/tau_F^2)) *asin(exp(-x^2/tau_F^2)))) *log(2 *u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)))-1)-(pi *log(u) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi) *((8 *x *exp(-x^2/tau_F^2))/(pi *tau_F^2 *sqrt(1-exp(-(2 *x^2)/tau_F^2)) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))+(16 *x *exp(-x^2/tau_F^2) *asin(exp(-x^2/tau_F^2)))/(pi^2 *tau_F^2 *sqrt(1-exp(-(2 *x^2)/tau_F^2)) *(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)^2)) *u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi*(1-(2 *asin(exp(-x^2/tau_F^2)))/pi))))/(2 *asin(exp(-x^2/tau_F^2)) *(2 *u^(-(4 *asin(exp(-x^2/tau_F^2)))/(pi*(1-(2 *asin(exp(-x^2/tau_F^2)))/pi)))-1)));
