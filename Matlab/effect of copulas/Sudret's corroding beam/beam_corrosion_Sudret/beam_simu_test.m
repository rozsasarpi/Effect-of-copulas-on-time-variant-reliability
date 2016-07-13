clear all
close all
clc

%time instant
%t = 10;
% correlation 'length' for F process
tau_F = 10/365;
% time increment
delta_t = tau_F*0.01;
% Pearsson correlation between 'adjactent' F realizations
ro_F = exp(-(delta_t/tau_F)^2);
% k_tau = exp(-(delta_t/tau_F));
% ro_F = sin(k_tau*pi/2);
%ro_F = 0;
% Kendall tau
k_tau = 2*asin(ro_F)/pi;

MU      = [10, 10];
SIGMA   = 1^2*ones(2,2);
SIGMA(~eye(2)) = SIGMA(~eye(2))*ro_F;
S       = 13;

copulatype = 'Gaussian';
% copulatype = 'Clayton';
% copulatype = 't';
% copulatype = 'Gumbel';
% copulatype = 'Frank';

%% direct integration solution
% marginal cdfs and pdfs
Fx1 = @(x) normcdf(x, MU(1), sqrt(SIGMA(1,1)));
Fx2 = @(x) normcdf(x, MU(2), sqrt(SIGMA(2,2)));
fx1 = @(x) normpdf(x, MU(1), sqrt(SIGMA(1,1)));
fx2 = @(x) normpdf(x, MU(2), sqrt(SIGMA(2,2)));
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
        % Placket copula
        %fx1x2 = @(x1, x2) theta*(1+(theta+1)*(Fx1(x1)+Fx2(x2)-2*Fx1(x1)*Fx2(x2)))/((1+(theta-1)*(Fx1(x1)+Fx2(x2)))^2-4*Fx1(x1)*Fx2(x2)*theta*(theta-1))^(3/2);
        fx1x2 = @(x1, x2) -theta.*(exp(-theta)-1).*exp(-theta.*(Fx1(x1)+Fx2(x2)))./((exp(-theta)-1)+(exp(-theta*Fx1(x1))-1).*(exp(-theta.*Fx2(x2))-1)).^2.*fx1(x1).*fx2(x2);
        
    case 'Gumbel'
        theta = 1/(1-k_tau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Gumbel', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
end

tic
Pf_d = integral2(@(x1,x2) fx1x2(x1,x2), S, 18, 5, S);
toc

% xx = 5:0.1:15;
% yy = xx;
% [X, Y] = meshgrid(xx,yy);
% Z = fx1x2(X,Y);
% surf(X,Y,Z)

% Rice formula for Gaussian processes, - normal marginals and multivariate normal joint distribution -> Gauss copula
disp('Rice''s formula:')
nu_p = 1/sqrt(2*pi)*sqrt(2*1/tau_F^2*theta)*normpdf(S, MU(1), sqrt(SIGMA(1,1)))
%Pf_R = nu_p*delta_t

%% simulate the random variables
tic
% number of simulations
N = 1e7;

x12     = mvnrnd(MU, SIGMA, N);

toc

% prob. that it is in the safe domain
id1 = x12(:,1) < S;
%std(id)/mean(id)/sqrt(N)

% prob. that it is in the failure domain
id2 = x12(:,2) > S;
ids = id1.*id2;
ns = sum(ids);
Pf_s = ns/N;

%nu_p = Pf_s/delta_t

disp('probability of failure:')
[Pf_d, Pf_s]
disp('outcrossing rate:')
[Pf_d, Pf_s]/delta_t