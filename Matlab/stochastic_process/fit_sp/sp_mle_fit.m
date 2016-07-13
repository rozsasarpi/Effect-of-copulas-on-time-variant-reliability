%Maximum likelihood fit of stochastic process with independent ensemble samples
% with normal marginal distribution!
%
%

function [parmhat, nll] = sp_mle_fit(X, sample, Options)

copula      = Options.Corr.copula;
likelihood  = Options.likelihood;

Corr_       = Options.Corr;

n_obs       = size(sample,1);
N           = size(sample,2);

x0          = [mean(sample(:)), std(sample(:)), 10];
lb          = [-Inf, 1e-3, 1e-3];
ub          = [Inf, Inf, Inf];

options     = optimoptions('fmincon');
options.Display = 'iter';

tic
switch lower(copula)
    case {'g', 'gauss'}
        [parmhat, nll] = fmincon(@nll_gauss,x0,[],[],[],[],lb,ub,[],options);
    case {'t', 'student'}
        [parmhat, nll] = fmincon(@nll_t,x0,[],[],[],[],lb,ub,[],options);
    case {'tev'}
%         [parmhat, nll] = fmincon(@nll_t,x0,[],[],[],[],lb,ub,[],options);
    case {'hr', 'husler-reiss', 'hüsler-reiss'}
        % vectors to interpolate Kendall tau and delta parameter
        load('hr_tau_delta.mat', 'hr_tau', 'hr_delta')
        [parmhat, nll] = fmincon(@nll_hr,x0,[],[],[],[],lb,ub,[],options);
    otherwise
        error(['Unknown copula function /copula/: ', copula])
end
toc

%--------------------------------------------------------------------------
% NEGATIVE LOGLIKELIHOOD FUNCTIONS
%--------------------------------------------------------------------------
%..........................................................................
% GAUSS COPULA
%..........................................................................
    function nll = nll_gauss(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        rho     = element_cov_matrix(X, Corr_);
        
        switch lower(likelihood)
            case {'f', 'full'}
                MU_     = ones(n_obs,1)*x(1);
                SIGMA_  = rho*x(2)^2;
                
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    nll     = nll + -log(mvnpdf(sample(:,ii), MU_,SIGMA_));
                    %         u       = normcdf(r(:,ii),x(1),x(2))';
                    %         nll     = nll + -log(copulapdf('Gaussian', u, rho)*prod(normpdf(r(:,ii),repmat(x(1),n_elem,1),repmat(x(2),n_elem,1))));
                end
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                C       = nchoosek(1:n_obs,2);
                M       = size(C,1);
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    % loop over the pairs in a particular ensemble
                    for jj = 1:M
                        MU_      = [x(1); x(1)];
                        r        = rho(C(jj,1),C(jj,2));
                        SIGMA_   = [1, r; r, 1]*x(2)^2;
                        %                        nll      = nll + -log(mvnpdf(sample(C(jj,:),ii), MU_, SIGMA_));
                        nll      = nll + -log(binormpdf(sample(C(jj,:),ii), MU_, SIGMA_));
                    end
                end
        end
    end

%..........................................................................
% T COPULA - WARNING! FIXED DOF!!
%..........................................................................
    function nll = nll_t(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        rho     = element_cov_matrix(X, Corr_);
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    u       = normcdf(sample(:,ii),x(1),x(2))';
                    nll     = nll + -log(copulapdf('t', u, rho, 2)*prod(normpdf(sample(:,ii),repmat(x(1),n_obs,1),repmat(x(2),n_obs,1))));
                end
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                C       = nchoosek(1:n_obs,2);
                M       = size(C,1);
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    % loop over the pairs in a particular ensemble
                    for jj = 1:M
                        r       = rho(C(jj,1),C(jj,2));
                        u       = normcdf(sample(C(jj,:),ii),x(1),x(2))';
%                         nll     = nll + -log(copulapdf('t', u, r, 2)*prod(normpdf(sample(C(jj,:),ii),repmat(x(1),2,1),repmat(x(2),2,1))));
                        nll     = nll + -log(bit_copulapdf(u, r, 2)*prod(normpdf(sample(C(jj,:),ii),repmat(x(1),2,1),repmat(x(2),2,1))));
                    end
                end
        end
    end

%..........................................................................
% HÜSLER-REISS COPULA
%..........................................................................
    function nll = nll_hr(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        rho     = element_cov_matrix(X, Corr_);
        % assuming that the covariance structure gives Kendall tau
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                error('Only pairwise likelihood is available for Hüsler-Reiss copula!')
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                C       = nchoosek(1:n_obs,2);
                M       = size(C,1);
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    % loop over the pairs in a particular ensemble
                    for jj = 1:M
                        r       = rho(C(jj,1),C(jj,2));
                        delta   = interp1(hr_tau, hr_delta, r);
                        u       = normcdf(sample(C(jj,:),ii), x(1), x(2))';
                        nll     = nll + -log(bihr_copulapdf(u, delta)*prod(normpdf(sample(C(jj,:),ii),repmat(x(1),2,1),repmat(x(2),2,1))));
                    end
                end
        end
    end

end