% Time-variant structural reliability
%
% g(t) = R - (G + S(t))
%
% stochastic process S(t) is inferred from actual snow measurements
% with various copulas (and marginals)
%
% RG = R-G with convolution => three-dimensional integral
%
%
% PROBLEM WITH TEV, tau to rho?
clearvars
close all
clc

copulas = {'gauss'};
marginals = {'gauss'};
ii = 1;
%--------------------------------------------------------------------------
% PROBABILISTIC MODELS
%--------------------------------------------------------------------------
Probvar.R.mean  = 250;
Probvar.R.cov   = 0.1;
Probvar.R.dist  = 2; %lognormal

Probvar.G.mean  = 60;
Probvar.G.cov   = 0.07;
Probvar.G.dist  = 1; %normal

% ML estimates, naive_fit.m
load(['ML_fit_', Corr.type, '_copn6_margn3.mat'], 'NLL', 'PARMHAT')
% par = [mean, std, corr_length]


Corr.length = par(3);
Corr.type   = 'gauss';
delta_t     = 0.01*par(3);

marginal    = marginals{ii};
copula      = copulas{ii};
%--------------------------------------------------------------------------
% CONVOLUTION
%--------------------------------------------------------------------------
[Probvar, pdRG] = conv_RG(Probvar);

frg         = @(rg) pdRG.fx_fun(rg);
% xx = linspace(pdRG.min, pdRG.max, 1e3);
% yy = pdRG.fx_fun(xx);
% plot(xx,yy)

% switch lower(marginal)
%     case {'g', 'gauss', 'norm', 'normal'}
%         marg_pdf = @(x, par) normpdf(x, par(:,1), par(:,2));
%         marg_cdf = @(x, par) normcdf(x, par(:,1), par(:,2));
%         
%     case {'lnorm', 'logn', 'lognorm', 'lognormal'}
%         marg_pdf = @(x, par) lognormpdf(x, par(:,1), par(:,2)./par(:,1));
%         marg_cdf = @(x, par) lognormcdf(x, par(:,1), par(:,2)./par(:,1));
%         
%     case {'gum', 'gumb', 'gumbel'}
%         marg_pdf = @(x, par) gumbelpdf(x, par(:,1), par(:,2));
%         marg_cdf = @(x, par) gumbelcdf(x, par(:,1), par(:,2)./par(:,1));
%     otherwise
%         error(['Unknown marginal distribution /marginal/: ', marginal])
% end

% stochastic load
Fx1 = @(x) marg_cdf(x, par);
Fx2 = @(x) marg_cdf(x, par);
fx1 = @(x) marg_pdf(x, par);
fx2 = @(x) marg_pdf(x, par);

CM          = cov_matrix([0, delta_t], Corr);
ktau        = CM(1,2);

switch lower(copula)
    case {'g', 'gauss'}
        rho   = sin(pi/2*ktau);
        theta = rho;
        
%         fx1x2 = @(x1, x2) reshape(binorm_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
        fx1x2 = @(x1, x2) reshape(copulapdf('Gauss',[reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
    case {'t', 'student'}
        rho   = sin(pi/2*ktau);
        theta = rho;
        
        fx1x2 = @(x1, x2) reshape(bit_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta, 2), size(x1,1), []).*fx1(x1).*fx2(x2);        
    case {'gum', 'gumb', 'gumbel'}
        theta   = 1/(1-ktau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Gumbel',[reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
        % Note: bigumb_copulapdf failes with 0/0 = NaN for cases e.g. high correlation
        % fx1x2 = @(x1, x2) reshape(bigumb_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
    case {'tev'}
        rho   = tev_ktau2rho(ktau);
        theta = rho;
        fx1x2 = @(x1, x2) reshape(bitev_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
    case {'hr'}
        delta = hr_ktau2delta(ktau);
        theta = delta;
        fx1x2 = @(x1, x2) reshape(bihr_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
    otherwise
        error(['Unknown copula function /copula/: ', copula])
end

% [X, Y] = meshgrid(-20:80, -20:80);
% Z = fx1x2(X,Y);
% surf(X,Y,Z)
% integral2(fx1x2, -50,300, -50,300)
%--------------------------------------------------------------------------
% ANALYSIS/CALCULATION
%--------------------------------------------------------------------------
% PF_0
fx1min = @(rg) rg ;
Pf_0 = integral2(@(rg,x1) frg(rg).*fx1(x1), -50,1000, fx1min,1000, 'AbsTol',1e-16, 'RelTol',1e-5);
disp(['Pf_0: ', num2str(Pf_0)])

% OUTCROSSING RATE
% integration limits
fx1max    = @(rg) rg;
fx2min    = @(rg,x1) rg;
frgx1x2   = @(rg,x1,x2) frg(rg).*fx1x2(x1,x2);

tic
% Pf_t = integral3(frgx1x2, -50,400 ,-50,fx1max, fx2min,400, 'AbsTol',1e-16, 'RelTol',1e-5);
Pf_t = integral3(frgx1x2, 0,400 ,400,fx1max, fx2min,0, 'AbsTol',1e-16, 'RelTol',1e-5); %WARNING
tint = toc;
disp(['Numerical integration:   ', num2str(tint), ' seconds.']);
disp('=====================')
Pf_t
-norminv(Pf_t)
nu_p = Pf_t/delta_t;
disp('Outcrossing rate, nu_p:')
nu_p

T = 50*365.25*0.2;
Pf = Pf_0 + T*nu_p
-norminv(Pf)