% Simple structure subjected to snow and wind (the other loads with resistance are 'integrated' into R' random variable)
% snow and wind are modelled as bivariate distribution
%
% SWE(X1)   - annual max snow water equivalent [mm]
% WS10(X2)  - accompanying (same day as annual swe max) horizontal wind speed measured 3 times a day [m/s]
% based on actual data from Budapest [E18.9, N47.5] 1961-2010
%
% limit state function:
% g = R' - (a_s*S + a_w*W)
%
% a_s and a_w are transformation parameters
%
% NOTE(S):
% parameters of the marginals and copulas (mle estiamtion) are given manually
% integration limit might need to be adjusted

clear all
% close all
clc

%==========================
% OPTIONS/ PARAMETERS
%==========================

% SWE, annual max, [mm], maximum likelihood fit
sigma_SWE   = 26.5748351;
mu_SWE      = 37.431650;
% transforamtion factor [N]
as          = 10;
% internal force from snow [kNm]
sigma_SWEas = as*sigma_SWE;
mu_SWEas    = as*mu_SWE;

% WS10, accompanying, [m/s], maximum likelihood fit
sigma_WS10  = 0.8182692;
mu_WS10     = 2.022557;
% transforamtion factor [N]
aw          = 200;
% internal force from snow [kNm]
sigma_WS10aw = aw*sigma_WS10;
mu_WS10aw    = aw*mu_WS10;

% copula
copulatype = 'Gaussian';
% copulatype = 't';
% copulatype = 'Clayton';
% copulatype = 'Gumbel';

% % copulatype = 'Frank';

% ro_F = 0.9;
% k_tau = 2*asin(ro_F)/pi;

%% direct integration solution
% marginal cdfs and pdfs
Fx1 = @(x) gevcdf(x, 0, sigma_SWEas, mu_SWEas);
Fx2 = @(x) gevcdf(x, 0, sigma_WS10aw, mu_WS10aw);
fx1 = @(x) gevpdf(x, 0, sigma_SWEas, mu_SWEas);
fx2 = @(x) gevpdf(x, 0, sigma_WS10aw, mu_WS10aw);

% resistance
fr = @(r) lognormpdf(r, 5000, 0.2); %kNm

% copula parameters are from maximum likelihood fit!
% WARNING! manually!
switch copulatype
    
    case 'Gaussian'
        theta = 0.3076610;%ro_F;
        fx1x2 = @(x1, x2) 1/sqrt(1-theta^2)*exp(-(norminv(Fx1(x1)).^2*theta^2 - 2*theta*norminv(Fx1(x1)).*norminv(Fx2(x2)) + norminv(Fx2(x2)).^2*theta^2)./(2*(1-theta^2))).*fx1(x1).*fx2(x2);
        %cpdf = @(x1, x2) reshape(copulapdf('Gaussian', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2); %for verification
    
    case 't' % dof = 2
        theta = 0.1642704;
        fx1x2 = @(x1, x2) reshape(copulapdf('t', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta, 2), size(x1,1), []).*fx1(x1).*fx2(x2);
    
    case 'Clayton'
        theta = 0.4384336;%2*k_tau/(1-k_tau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Clayton', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
        
%     case 'Frank'
%         theta = fzero(@(alpha) copulastat('Frank', alpha, 'type', 'Kendall') - k_tau, ro_F);
%         %fx1x2 = @(x1, x2) reshape(copulapdf('Frank', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);     
%         fx1x2 = @(x1, x2) -theta.*(exp(-theta)-1).*exp(-theta.*(Fx1(x1)+Fx2(x2)))./((exp(-theta)-1)+(exp(-theta*Fx1(x1))-1).*(exp(-theta.*Fx2(x2))-1)).^2.*fx1(x1).*fx2(x2);
%     case 'Plackett'
%         % not yet tested
%         %fx1x2 = @(x1, x2) theta*(1+(theta+1)*(Fx1(x1)+Fx2(x2)-2*Fx1(x1)*Fx2(x2)))/((1+(theta-1)*(Fx1(x1)+Fx2(x2)))^2-4*Fx1(x1)*Fx2(x2)*theta*(theta-1))^(3/2);
    case 'Gumbel'
        theta = 1.1557412;
        fx1x2 = @(x1, x2) reshape(copulapdf('Gumbel', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
end

% % plot
% x = 0:9000;
% plot(x,fx1(x))
% hold on
% plot(x,fx2(x),'r')
% plot(x,fr(x),'g--')
% legend('SWE','WS10','R´')

%==========================
% ANALYSIS/CALCULATION
%==========================
% integration limits
fx2min    = @(r,x1) r - x1;
frx1x2    = @(r,x1,x2) fr(r).*fx1x2(x1,x2);
tic
% WARNING! integration limits!
Pf = integral3(frx1x2, 100,7000 ,100,7000, fx2min,6000, 'AbsTol',1e-12, 'RelTol',1e-5);
% Pf = integral3(frx1x2, -1000,7000 ,-1000,7000, fx2min,7000, 'AbsTol',1e-12, 'RelTol',1e-5);
% Pf = integral3(frx1x2, 2,6 ,-0.01,2, fx2min,7, 'AbsTol',1e-12, 'RelTol',1e-5);
tint = toc;
disp(['Numerical integration:   ', num2str(tint), ' seconds.']);
disp('=====================')

% Pf_d

disp(['Probability of failure: ', num2str(Pf)])
disp(['Reliability index: ', num2str(-norminv(Pf))])