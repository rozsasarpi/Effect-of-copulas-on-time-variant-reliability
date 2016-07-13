% Fit copula to snow-wind multivariate extremes
% comparison to R results

clear all
close all
clc

% results using R (some checked with two different built-in functions)
%              theta        AIC         BIC
% Gauss    0.3076610 -1.8376807  0.07434229
% t        0.1642704  6.7935750  8.70559796
% Clayton  0.4384336 -2.2949996 -0.38297658
% Gumbel   1.1557412  0.3377976  2.24982059
% Plackett 1.8311975 -0.3479257  1.56409731
% Frank    1.3936066 -0.6547995  1.25722348

% Gauss, t, Clayton and Gumbel theta parameters are the same as those used in matlab, thus directly comparable

load('annual_snow_max_acc_wind.mat')

% to test on other data
% rng(111)
% R = mvnrnd([0,0], [1,0.5;0.5,1],10);
% 
% snow = R(:,1);
% wind = R(:,2);
% 
% ro = corr(snow, wind);

% pseudo observations
u_snow = tiedrank(snow)/(length(snow)+1);
u_wind = tiedrank(wind)/(length(wind)+1);

plot(u_snow, u_wind, 'o')
xlabel('u_{snow}')
ylabel('u_{wind}')
axis([0,1,0,1])
% axis equal

% ====================================================
% Copula MLE fit
% ====================================================
disp('Maximum likelihood copula fit.')
disp('bi - Matlab built-in')
disp('hm - home-made, direct optimization')
fprintf('\n');fprintf('\n')

% ....................................................
% Gaussian
% ....................................................
disp('Gaussian copula')
% built-in function - ! it looks like that the copulafit is not maximum likelihood approach for gaussian copula..
% something between Kendall' tau and Spearman's rho..
%options     = statset('TolX', 1e-7, 'MaxIter', 1000);
% options 'is not applicable to the ''Gaussian'' family'! ridiculous..
cov_mx      = copulafit('Gaussian', [u_snow, u_wind]);
ro_gauss_bi = cov_mx(1,2);

% neagative log-likelihood
nLL = @(x) -sum(log(copulapdf('Gaussian', [u_snow, u_wind], x)));

% home-made solu
ro_gauss_hm = fminbnd(nLL, -1, 1);

tb_gauss = table(ro_gauss_bi, ro_gauss_hm);
disp(tb_gauss)

% ....................................................
% t copula, df = 2,fixed
% ....................................................
disp('t copula')
% no built-in function for fixed df

% neagative log-likelihood
nLL = @(x) -sum(log(copulapdf('t', [u_snow, u_wind], x, 2)));

% home-made solu
ro_t_hm = fminbnd(nLL, -1, 1);

tb_t = table(ro_t_hm);
disp(tb_t)
