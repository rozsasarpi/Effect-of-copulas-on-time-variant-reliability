% Maximum likelihood fit of stochastic process


function naive_fit

close all
clc

% rng(333) % for reproducibility
% copulas = {'gauss'};
% marginals = {'lognormal'};
copulas     = {'gauss', 't', 'gumbel', 'rotgumbel', 'rotclayton'};
marginals   = {'gauss', 'lognormal', 'gumbel'};
n           = length(copulas);
m           = length(marginals);
NLL         = nan(n,m);
PARMHAT     = nan(n,m,3);
for ii = 1:n
    for jj = 1:m
        %--------------------------------------------------------------------------
        % OPTIONS
        %--------------------------------------------------------------------------
        
        % N           = 50;           % number of ensembles
        % X           = 1:10;         % distance of observations within an ensemble
        Corr.length = 10/365;           % correlation length (used in the autocorrelation function)
        Corr.pow    = 2;            % parameter used in the autocorrelation function
        Corr.type   = 'cauchy';
        Corr.copula = copulas{ii};
        likelihood  = 'p';
        marginal    = marginals{jj};
        
        % mu          = 10;           % mean of the marginal of intensity
        % sigma       = 2;            % std of the marginal of intensity
        
        Options.Corr        = Corr;
        Options.likelihood  = likelihood;
        Options.marginal    = marginal;
        %--------------------------------------------------------------------------
        % GENERATE SYNTHETIC OBSERVATIONS / RANDOM SAMPLES
        %--------------------------------------------------------------------------
        % n_obs       = length(X);
        % cov_matrix  = element_cov_matrix(X, Corr);
        % sample      = mvnrnd(ones(n_obs,1)*mu, cov_matrix*sigma^2, N)';
        % sample      = sample + normrnd(0,1,length(X),N); %noise
        %
        % sample      = mat2cell(sample,length(X), ones(1,N,1))';
        % X           = repmat(mat2cell(X', size(X,2), size(X,1)),N,1);
        load('D:\Working folder\Matlab working folder\time variant reliability\phi2 method\effect of copulas\stochastic_snow_model\snow_ensembles.mat')
        % convert time to years
        X = cellfun(@(x) x./365, X, 'UniformOutput', false);
        
        % parmhat = [mean, std, corr_length]
        % [parmhat, nll] = sp_mle_fit(X, sample, Options);
        [parmhat, nll] = sp_mle_fit_cell(X, sample, Options);
        
        NLL(ii,jj) = nll;
        PARMHAT(ii,jj,:) = parmhat;
        
    end
end
% nll
% parmhat
% parmhat(3)*365

% n = size(X);
% for ii = 1:n
%     if ii == 1
%         xx = X{ii};
%         ss = sample{ii};
%     else
%         xx = [xx; X{ii}];
%         ss = [ss; sample{ii}];
%     end
%     plot(X{ii}, sample{ii}, 'o-b')
%     hold on
% end

% plot(xx,ss)

% plot(sample)
% xt = min(xlim);
% yt = max(ylim) - 0.05*diff(ylim);
% text(xt,yt,['Copula: ', Corr.copula])
% yt = max(ylim) - 0.10*diff(ylim);
% text(xt,yt,['Autocorr.: ', Corr.type])
% yt = max(ylim) - 0.15*diff(ylim);
% text(xt,yt,['Likelihood.: ', Options.likelihood])
save(['ML_fit_', Corr.type, '_copn', num2str(n), '_margn', num2str(m), '.mat'], 'NLL', 'PARMHAT')
end