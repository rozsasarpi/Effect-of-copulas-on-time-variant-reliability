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

mp.Digits(100);

copulas         = {'gauss', 't', 'gumbel', 'rotgumbel', 'rotclayton'};
% copulas         = {'gauss', 't'};
% copulas         = {'t', 'gumbel', 'rotgumbel', 'rotclayton'};
marginals       = {'gauss', 'lognormal', 'gumbel'};
acorr_fun       = 'cauchy';

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
load(['D:\Working folder\Matlab working folder\stochastic_process\fit_sp\ML_fit_', acorr_fun, '_copn5_margn3.mat'], 'NLL', 'PARMHAT')
copula_pat  = {'gauss', 't', 'gumbel', 'rotgumbel', 'rotclayton'};
marginal_pat= {'gauss', 'lognormal', 'gumbel'};

n           = length(copulas);
m           = length(marginals);
PF_t        = nan(n,m);
NU_p        = nan(n,m);
BETA        = nan(n,m);
for ii = 1:n
    for jj = 1:m
        
        % par = [mean, std, corr_length]
        parmhat         = PARMHAT(1,1,:);
        
        marginal        = marginal_pat(strcmpi(marginal_pat, marginals{jj}));
        marginal        = marginal{:};
        copula          = copula_pat(strcmpi(copula_pat, copulas{ii}));
        copula          = copula{:};
        %--------------------------------------------------------------------------
        % CONVOLUTION, RG = R - G
        %--------------------------------------------------------------------------
        [Probvar, pdRG] = conv_RG(Probvar);
        
        xx = 0:0.1:500;
        yy = pdRG.fx_fun(xx);
        % plot(xx,yy)
        fr              = @(r) interp1(xx, yy, r);
        % fr              = @(r) pdRG.fx_fun(r);
        
        tau_F           = parmhat(3);
        
        % time increment
        delta_t         = tau_F*mp('0.001');
        % Pearson correlation between 'adjactent' S realizations
        % #########################################################################
        % Gaussian
        ro_F            = exp(-(delta_t/tau_F)^2);
        % Cauchy
        %     ro_F    = (1+(delta_t/tau_F)^2).^(-2);
        % #########################################################################
        % Kendall tau
        k_tau           = 2*asin(ro_F)/mp('pi');
        
        % stochastic load
        MU              = [parmhat(1), parmhat(1)];
        SIGMA           = parmhat(2)^2*ones(2,2);
        SIGMA(~eye(2))  = SIGMA(~eye(2))*ro_F;
        
        % integration limits
        int_l           = mp('50');
        int_u           = mp('200');
        
        
        % copula
        switch lower(copula)
            case {'gaussian', 'gauss'}
                theta = ro_F;
                F     = @(r) binorm_rect_mp(int_l, r, r, int_u, MU, SIGMA, marginal).*fr(r);
            case 't' % dof = 2
                %theta = ro_F;
                theta = sin(k_tau*pi/2);
                F     = @(r) bit_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta, 2, marginal).*fr(r);
            case 'rotclayton'
                theta = 2*k_tau/(1-k_tau);
                F     = @(r) birotclay_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta, marginal).*fr(r);
            case 'gumbel'
                theta = 1/(1-k_tau);
                F     = @(r) bigumb_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta, marginal).*fr(r);
            case 'rotgumbel'
                theta = 1/(1-k_tau);
                F     = @(r) birotgumb_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta, marginal).*fr(r);
            otherwise
                error(['Unknown copula type:', copulatype])
        end
        
        %==========================
        % ANALYSIS/CALCULATION
        %==========================
        tic
        Pf_t    = quadgk(F, int_l, int_u, 'RelTol',mp('1e-8'),'AbsTol',mp('1e-12'))
        toc
        Pf_0    = 0;
        nu_p    = Pf_t/delta_t
        
        T       = 1*365.25*0.2;
        Pf      = Pf_0 + T*nu_p;
        beta    = -norminv(double(Pf))
        
        PF_t(ii,jj) = Pf_t;
        NU_p(ii,jj) = nu_p;
        BETA(ii,jj) = beta;
        
    end
end

save(['reli_results_',acorr_fun,'.mat'], 'PF_t', 'NU_p', 'BETA')
