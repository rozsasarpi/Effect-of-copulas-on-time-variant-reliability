clear all
close all
clc

%--------------------------------------------------------------------------
tau_F = 1;
delta_t = (0:0.01:2)*tau_F;

ro_F = exp(-(delta_t/tau_F).^2);
k_tau = 2*asin(ro_F)/pi;

subplot(1,2,1)
plot(delta_t, ro_F)
hold on
plot(delta_t, k_tau)

%
% k_tau = exp(-(delta_t/tau_F).^2);
k_tau = 2/pi*asin(exp(-(delta_t/tau_F).^2));
ro_F = sin(k_tau*pi/2);

subplot(1,2,2)
plot(delta_t, ro_F)
hold on
plot(delta_t, k_tau)


%--------------------------------------------------------------------------
% Kendall tau and copula parameter (theta) connection
%--------------------------------------------------------------------------

k_tau = -1:0.01:1;
theta_ga    = 2/pi*asin(k_tau);
theta_t     = theta_ga;
theta_c     = 2*k_tau./(1-k_tau);
theta_gu    = 1./(1-k_tau);

figure
plot(k_tau, [theta_ga', theta_c', theta_gu'])

ylabel('\theta')
xlabel('\tau_K')
ylim([-2,10])
legend('gauss=t', 'clayton', 'gumbel', 'Location', 'Northwest')