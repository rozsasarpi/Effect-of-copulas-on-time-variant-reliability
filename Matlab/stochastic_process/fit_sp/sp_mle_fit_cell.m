%Maximum likelihood fit of stochastic process with independent ensemble samples
% with normal marginal distribution!
%
%

function [parmhat, nll] = sp_mle_fit_cell(X, sample, Options)

copula      = Options.Corr.copula;
marginal    = Options.marginal;
likelihood  = Options.likelihood;

Corr_       = Options.Corr;

n_obs       = cellfun(@(x) length(x), sample);
N           = size(sample,1);

x0          = [mean(cellfun(@mean, sample)), mean(cellfun(@std, sample)), 10/365];
lb          = [1e-3, 1e-3, 1e-3]; % constrain on mean to prevent lognormal jumping to negative region - symmetric
ub          = [Inf, Inf, Inf];

options     = optimoptions('fmincon');
options.Display = 'iter';
% options.TolX =

% MARGINAL DISTRIBUTION
switch lower(marginal)
    case {'g', 'gauss', 'norm', 'normal'}
        marg_pdf = @(x, par) normpdf(x, par(:,1), par(:,2));
        marg_cdf = @(x, par) normcdf(x, par(:,1), par(:,2));
        
    case {'lnorm', 'logn', 'lognorm', 'lognormal'}
        marg_pdf = @(x, par) lognormpdf(x, par(:,1), par(:,2)./par(:,1));
        marg_cdf = @(x, par) lognormcdf(x, par(:,1), par(:,2)./par(:,1));
        
    case {'gum', 'gumb', 'gumbel'}
        marg_pdf = @(x, par) gumbelpdf(x, par(:,1), par(:,2));
        marg_cdf = @(x, par) gumbelcdf(x, par(:,1), par(:,2)./par(:,1));
    otherwise
        error(['Unknown marginal distribution /marginal/: ', marginal])
end
tic
% COPULA FUNCTION
switch lower(copula)
    case {'g', 'gauss'}
        %         fun = @(x) nll_gauss([x0(1), x0(2), x]);
        %         [parmhat, nll] = fmincon(fun,x0(3),[],[],[],[],lb,ub,[],options);
        [parmhat, nll] = fmincon(@nll_gauss,x0,[],[],[],[],lb,ub,[],options);
    case {'t', 'student'}
        %         fun = @(x) nll_t([x0(1), x0(2), x]);
        %         [parmhat, nll] = fmincon(fun,x0(3),[],[],[],[],lb,ub,[],options);
        [parmhat, nll] = fmincon(@nll_t,x0,[],[],[],[],lb,ub,[],options);
    case {'gum', 'gumb', 'gumbel'}
        [parmhat, nll] = fmincon(@nll_gumb,x0,[],[],[],[],lb,ub,[],options);
    case {'rotgum', 'rotgumb', 'rotgumbel'}
        [parmhat, nll] = fmincon(@nll_rotgumb,x0,[],[],[],[],lb,ub,[],options);
    case {'rotclay', 'rotclayton'}
        [parmhat, nll] = fmincon(@nll_rotclay,x0,[],[],[],[],lb,ub,[],options);
    case {'tev'}
        [parmhat, nll] = fmincon(@nll_tev,x0,[],[],[],[],lb,ub,[],options);
    case {'hr', 'husler-reiss', 'hüsler-reiss'}
        [parmhat, nll] = fmincon(@nll_hr,x0,[],[],[],[],lb,ub,[],options);
    otherwise
        error(['Unknown copula function /copula/: ', copula])
end
toc
% x1 = -25:2:25;
% nn = length(x1);
% y1 = nan(nn,1);
% for kk = 1:nn
%     y1(kk) = nll_gauss([x1(kk), parmhat(2), parmhat(3)]);
% end
% plot(x1, y1)
% hold on
% plot(parmhat(1), nll, 'o')
% 
% x1 = 15:2:45;
% nn = length(x1);
% y1 = nan(nn,1);
% for kk = 1:nn
%     y1(kk) = nll_gauss([parmhat(1), x1(kk), parmhat(3)]);
% end
% figure
% plot(x1, y1)
% hold on
% plot(parmhat(2), nll, 'o')
% xlabel('parmhat(2)')
% ylabel('-log(L)')
%--------------------------------------------------------------------------
% NEGATIVE LOGLIKELIHOOD FUNCTIONS
%--------------------------------------------------------------------------
%..........................................................................
% GAUSS COPULA
%..........................................................................
    function nll = nll_gauss(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        
        switch lower(likelihood)
            case {'f', 'full'}
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    ktau    = cov_matrix(X{ii}, Corr_);
                    rho     = sin(ktau*pi/2);
                    u       = marg_cdf(sample{ii},[x(1),x(2)])';
                    nll     = nll + -log(copulapdf('Gaussian', u, rho)*prod(marg_pdf(sample{ii},[repmat(x(1),n_obs(ii),1),repmat(x(2),n_obs(ii),1)])));
                end
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            rho     = sin(ktau*pi/2);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            nll     = nll + -log(binorm_copulapdf(u, rho)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
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
        
        switch lower(likelihood)
            case {'f', 'full'}
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    ktau    = cov_matrix(X{ii}, Corr_);
                    rho     = sin(ktau*pi/2);
                    u       = marg_cdf(sample{ii},[x(1),x(2)])';
                    nll     = nll + -log(copulapdf('t', u, rho, 2)*prod(marg_pdf(sample{ii},[repmat(x(1),n_obs(ii),1),repmat(x(2),n_obs(ii),1)])));
                end
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            rho     = sin(ktau*pi/2);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            %                         nll     = nll + -log(copulapdf('t', u, r, 2)*prod(normpdf(sample(C(jj,:),ii),repmat(x(1),2,1),repmat(x(2),2,1))));
                            nll     = nll + -log(bit_copulapdf(u, rho, 2)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
                    end
                end
        end
    end

%..........................................................................
% GUMBEL COPULA
%..........................................................................
    function nll = nll_gumb(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                error('Only pairwise likelihood is available for Gumbel copula!')
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        % assuming that the covariance structure gives Kendall tau
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            % assuming that the covariance structure gives Kendall tau
                            theta   = 1/(1-ktau);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            nll     = nll + -log(bigumb_copulapdf(u, theta)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
                    end
                end
        end
    end

%..........................................................................
% ROTATED (180°) GUMBEL COPULA
%..........................................................................
    function nll = nll_rotgumb(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                error('Only pairwise likelihood is available for rotated Gumbel copula!')
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        % assuming that the covariance structure gives Kendall tau
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            % assuming that the covariance structure gives Kendall tau
                            theta   = 1/(1-ktau);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            nll     = nll + -log(bigumb_copulapdf(1-u, theta)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
                    end
                end
        end
    end

%..........................................................................
% ROTATED (180°) CLAYTON COPULA
%..........................................................................
    function nll = nll_rotclay(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                error('Only pairwise likelihood is available for rotated Gumbel copula!')
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        % assuming that the covariance structure gives Kendall tau
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            % assuming that the covariance structure gives Kendall tau
                            theta   = 1/(1-ktau);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            nll     = nll + -log(biclay_copulapdf(1-u, theta)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
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
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                error('Only pairwise likelihood is available for Hüsler-Reiss copula!')
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        % assuming that the covariance structure gives Kendall tau
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            % assuming that the covariance structure gives Kendall tau
                            delta   = hr_ktau2delta(ktau);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            nll     = nll + -log(bihr_copulapdf(u, delta)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
                    end
                end
        end
    end

%..........................................................................
% EXTREME T COPULA - WARNING! FIXED DOF!!
%..........................................................................
    function nll = nll_tev(x)
        % x = [mean, std, corr_length]
        Corr_.length = x(3);
        
        switch lower(likelihood)
            case {'f', 'full'}
                
                error('Only pairwise likelihood is available for tev copula!')
                
            case {'p', 'pairwise'}  %quasi-maximum likelihood
                nll     = 0;
                % loop over the ensembles
                for ii = 1:N
                    s       = sample{ii};
                    if n_obs(ii) == 1
                        nll =  nll + -log(marg_pdf(s, [x(1), x(2)]));
                    else
                        % assuming that the covariance structure gives Kendall tau
                        Ktau    = cov_matrix(X{ii}, Corr_);
                        C       = nchoosek(1:n_obs(ii),2);
                        M       = size(C,1);
                        % loop over the pairs in a particular ensemble
                        for jj = 1:M
                            ktau    = Ktau(C(jj,1),C(jj,2));
                            rho     = tev_ktau2rho(ktau);
                            u       = marg_cdf(s(C(jj,:)),[x(1),x(2)])';
                            nll     = nll + -log(bitev_copulapdf(u, rho, 2)*prod(marg_pdf(s(C(jj,:)),[repmat(x(1),2,1),repmat(x(2),2,1)])));
                        end
                    end
                end
        end
    end

end