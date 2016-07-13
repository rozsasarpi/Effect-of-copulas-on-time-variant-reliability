% Time-variant reliability analysis with continuous stochastic process
% investigation of the effect of copula assumption
%
% g = R - S(t)
%
% R    - constant
% S(t) - continuous stochastic process

clearvars
close all
clc

%--------------------------------------------------------------------------
% OPTIONS
%--------------------------------------------------------------------------

% COPULA
% copulatype = {'Gaussian', 't', 'Clayton', 'Gumbel'};
copulatype = {'Gaussian'};
% copulatype = {'t'};
% copulatype = {'Clayton'};
% copulatype = {'Gumbel'};

% copulatype = 'Frank'; % NOT WORKING

% AUTOCORRELATION
Corr.type = 'g';
Corr.pow  = 2;

% RESISTANCE
% resistance, to reach the desired Pf0 value, has no effect, just for
% convinience
R = 354.45;%339.91;%290.135;%640.7;

% correlation 'length' for S process
tau_Fv = 1/365;%[0.1, 1, 10, 100, 1000]/365;

%--------------------------------------------------------------------------
% ANALYSIS/CALCULATION
%--------------------------------------------------------------------------
mm  = numel(copulatype);
nn  = numel(tau_Fv);

Pf  = nan(mm, nn);
Pf2 = nan(mm, nn);
% loop over copula types
for jj = 1:mm
    
    copula = copulatype{jj};
    
    % loop over correlation lengths
    for ii = 1:nn
        
        tau_F   = tau_Fv(ii);
        Corr.length = tau_F;
        % time increment
        frac    = 0.1;
        delta_t = tau_F*frac;
        
        % ktau correlation between 'adjactent' F realizations
%         k_tau   = 2/pi*asin(exp(-(delta_t/tau_F).^2)); % to get back the typical correlation function for pearson rho
                CM = cov_matrix([0, delta_t], Corr);
                k_tau = 2/pi*asin(CM(1,2));
        
        % Pearson correlation for Gauss and t copulas
        ro_F    = sin(k_tau*pi/2);
        
        % marginal properties, have no effect, just for convinience
        meanS   = 100;
        covS    = 0.4;
        MU      = [meanS, meanS];
        SIGMA   = (covS*meanS)^2*ones(2,2);
        SIGMA(~eye(2)) = SIGMA(~eye(2))*ro_F;
        
        %------------------------------------------------------------------
        % DIRECT INTEGRATION - INITIALIZATION
        %------------------------------------------------------------------
        % marginal cdfs and pdfs
        % stochastic load
        Fx1 = @(x) normcdf(x, MU(1), sqrt(SIGMA(1,1)));
        Fx2 = @(x) normcdf(x, MU(2), sqrt(SIGMA(2,2)));
        fx1 = @(x) normpdf(x, MU(1), sqrt(SIGMA(1,1)));
        fx2 = @(x) normpdf(x, MU(2), sqrt(SIGMA(2,2)));
        
        % Fx1 = @(x) lognormcdf(x, MU(1), sqrt(SIGMA(1,1))/MU(1));
        % Fx2 = @(x) lognormcdf(x, MU(2), sqrt(SIGMA(2,2))/MU(2));
        % fx1 = @(x) lognormpdf(x, MU(1), sqrt(SIGMA(1,1))/MU(1));
        % fx2 = @(x) lognormpdf(x, MU(2), sqrt(SIGMA(2,2))/MU(2));
        
        
        % copula
        switch copula
            
            case 'Gaussian'
                theta = ro_F;
%                 fx1x2 = @(x1, x2) 1/sqrt(1-theta^2)*exp(-(norminv(Fx1(x1)).^2*theta^2 - 2*theta*norminv(Fx1(x1)).*norminv(Fx2(x2)) + norminv(Fx2(x2)).^2*theta^2)./(2*(1-theta^2))).*fx1(x1).*fx2(x2);
                fx1x2 = @(x1, x2) reshape(binorm_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
                %fx1x2 = @(x1, x2) reshape(copulapdf('Gaussian', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2); %for verification
                
            case 't' % dof = 2
                theta = ro_F;
                
%                 fx1x2 = @(x1, x2) reshape(copulapdf('t', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta, 2), size(x1,1), []).*fx1(x1).*fx2(x2);
                fx1x2 = @(x1, x2) reshape(bit_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta, 2), size(x1,1), []).*fx1(x1).*fx2(x2);
            case 'Clayton'
                theta = 2*k_tau/(1-k_tau);
                fx1x2 = @(x1, x2) reshape(copulapdf('Clayton', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
%                 fx1x2 = @(x1, x2) reshape(biclay_copulapdf([reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
            case 'Frank'
                theta = fzero(@(alpha) copulastat('Frank', alpha, 'type', 'Kendall') - k_tau, ro_F);
                %fx1x2 = @(x1, x2) reshape(copulapdf('Frank', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
                fx1x2 = @(x1, x2) -theta.*(exp(-theta)-1).*exp(-theta.*(Fx1(x1)+Fx2(x2)))./((exp(-theta)-1)+(exp(-theta*Fx1(x1))-1).*(exp(-theta.*Fx2(x2))-1)).^2.*fx1(x1).*fx2(x2);
            case 'Plackett'
                % not yet tested
                %fx1x2 = @(x1, x2) theta*(1+(theta+1)*(Fx1(x1)+Fx2(x2)-2*Fx1(x1)*Fx2(x2)))/((1+(theta-1)*(Fx1(x1)+Fx2(x2)))^2-4*Fx1(x1)*Fx2(x2)*theta*(theta-1))^(3/2);
            case 'Gumbel'
                theta = 1/(1-k_tau);
                fx1x2 = @(x1, x2) reshape(copulapdf('Gumbel', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
        end
        
        %------------------------------------------------------------------
        % DIRECT INTEGRATION - CALCULATION - BUILT-IN COPULACDF
        %------------------------------------------------------------------
        % survival probability of a component
        Ps_comp = normcdf(R, meanS, covS*meanS);
        if strcmp(copula, 't')
            Pf_sys  = Ps_comp - copulacdf(copula, [Ps_comp, Ps_comp], theta, 2);
        else
            Pf_sys  = Ps_comp - copulacdf(copula, [Ps_comp, Ps_comp], theta);
        end
        
        %------------------------------------------------------------------
        % DIRECT INTEGRATION - CALCULATION - BUILT-IN INTEGRAL
        %------------------------------------------------------------------
        % numerical integration check
        % integration limits
        tic
        Pf_sys2 = integral2(fx1x2, 0,R, R,400, 'AbsTol',1e-16, 'RelTol',1e-5); %AbsTol !!
        tint = toc;
        disp(['Numerical integration:   ', num2str(tint), ' seconds.']);
        disp('=====================')
        
        nu_p = Pf_sys/delta_t;
        nu_p2 = Pf_sys2/delta_t;
        
        Pf0 = 1 - Ps_comp;
        Pf(jj,ii) = Pf0 + nu_p*50;
        Pf2(jj,ii) = Pf0 + nu_p2*50;
        
    end
end
copulatype
Pf0
Pf
Pf2
-norminv(Pf)

bsxfun(@rdivide, Pf, Pf(1,:))