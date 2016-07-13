% g = R - S(t)
%
% R    - constant
% S(t) - continuous stochastic process

clear all
close all
clc

%==========================
% OPTIONS
%==========================


% correlation 'length' for F process
tau_Fv = [0.1, 1, 10, 100, 1000]/365;

Pf = zeros(numel(tau_Fv),1);
for i = 1:numel(tau_Fv)
    tau_F = tau_Fv(i);
    % time increment
    delta_t = tau_F*0.1;
    % Pearson correlation between 'adjactent' F realizations
    ro_F = exp(-(delta_t/tau_F)^2);
    % k_tau = exp(-(delta_t/tau_F));
    % ro_F = sin(k_tau*pi/2);
%     ro_F = 0.1;%1-1e-16;
    % Kendall tau
    k_tau = 2*asin(ro_F)/pi;
    

    meanS   = 100;
    covS    = 0.4;
    MU      = [meanS, meanS];
    SIGMA   = (covS*meanS)^2*ones(2,2);
    SIGMA(~eye(2)) = SIGMA(~eye(2))*ro_F;
    
%     copulatype = 'Gaussian';
%     copulatype = 't';
    copulatype = 'Clayton';
%     copulatype = 'Gumbel';
    
    % copulatype = 'Frank';
    
    %% direct integration solution
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
    
    % resistance
    R =339.91;%290.135;%640.7;

    % copula
    switch copulatype
        
        case 'Gaussian'
            theta = ro_F;
            %         theta = exp(-(delta_t/tau_F)^2);%sin(k_tau*pi/2);
            fx1x2 = @(x1, x2) 1/sqrt(1-theta^2)*exp(-(norminv(Fx1(x1)).^2*theta^2 - 2*theta*norminv(Fx1(x1)).*norminv(Fx2(x2)) + norminv(Fx2(x2)).^2*theta^2)./(2*(1-theta^2))).*fx1(x1).*fx2(x2);
            %cpdf = @(x1, x2) reshape(copulapdf('Gaussian', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2); %for verification
            
        case 't' % dof = 2
            %theta = ro_F;
            theta = sin(k_tau*pi/2);
            %cpdf = @(x1, x2)
            fx1x2 = @(x1, x2) reshape(copulapdf('t', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta, 2), size(x1,1), []).*fx1(x1).*fx2(x2);
            
        case 'Clayton'
            theta = 2*k_tau/(1-k_tau);
            fx1x2 = @(x1, x2) reshape(copulapdf('Clayton', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
            
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
    
    %==========================
    % ANALYSIS/CALCULATION
    %==========================
    % survival probability of a component
    Ps_comp = normcdf(R,meanS,covS*meanS);
    if strcmp(copulatype, 't')
        Pf_sys  = Ps_comp-copulacdf(copulatype, [Ps_comp, Ps_comp], theta, 2);
    else
        Pf_sys  = Ps_comp-copulacdf(copulatype, [Ps_comp, Ps_comp], theta);
    end
    
    % numerical integration check
    % integration limits
    % tic
    % Pf_sys2 = integral2(fx1x2, 0,R, R,400, 'AbsTol',1e-16, 'RelTol',1e-5);
    % tint = toc;
    % disp(['Numerical integration:   ', num2str(tint), ' seconds.']);
    % disp('=====================')
    
%     tau_F*365
    % Pf_sys
    % Pf_sys2
    nu_p = Pf_sys/delta_t;
    
    % initial probability of failure, t=0
    % frx1 = @(r,x1) fr(r).*fx1(x1);
    % Pf_0 = integral2(frx1, 500,9500 ,fx1max,9500, 'AbsTol',1e-12, 'RelTol',1e-5);
    
    % disp('nu_p:')
%     nu_p'
    
    Pf0 = 1-Ps_comp;
    Pf(i) = Pf0 + nu_p*50;
%     Pf(i)
%     -norminv(Pf(i))
%     disp('=====================')
end
Pf0
Pf'
-norminv(Pf)'