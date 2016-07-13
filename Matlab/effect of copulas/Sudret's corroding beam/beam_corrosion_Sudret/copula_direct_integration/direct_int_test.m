% Sudret et al. (2002). Comparison of methods for computing the probability of failure in time-variant reliability using the outcrossing approach
% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
% simply supported beam subjected to stochastic load (F) and corrosion

clearvars
close all
clc

%==========================
% OPTIONS
%==========================
tt = [0, 2, 5, 10, 15, 20];
sim = 0;        % calculate the pdf of the R-Msw random variable?, 0 - no; 1 - yes
% tt = 0;
% frac_t = [0.1, 0.01, 0.001];
% frac_t = [0.001, 0.0001];
frac_t = 0.1;
mm = length(frac_t);
for jj = 1:mm
nu_p = zeros(numel(tt),1);
for i = 1:numel(tt)
%time instant
t = tt(i);
% correlation 'length' for F process
tau_F = 1/365;
% span
L       = 5;
% time increment
delta_t = tau_F*frac_t(jj);
% Pearsson correlation between 'adjactent' F realizations
ro_F = exp(-(delta_t/tau_F)^2);
% k_tau = exp(-(delta_t/tau_F));
% ro_F = sin(k_tau*pi/2);
%ro_F = 0;
% Kendall tau
k_tau = 2*asin(ro_F)/pi;

% ro_F = 0.5;
MU      = [3500, 3500]*L/4;
SIGMA   = (L/4*3500*0.2)^2*ones(2,2);
SIGMA(~eye(2)) = SIGMA(~eye(2))*ro_F;

copulatype = 'Gaussian';
% copulatype = 't';
% copulatype = 'rotClayton';% <0.01 step is required for convergence of outcrossing rate!
% copulatype = 'Gumbel';
% copulatype = 'rotGumbel'; % <0.001 step is required for convergence of outcrossing rate!

% copulatype = 'Frank';

%% direct integration solution
% marginal cdfs and pdfs
% stochastic load
Fx1 = @(x) normcdf(x, MU(1), sqrt(SIGMA(1,1)));
Fx2 = @(x) normcdf(x, MU(2), sqrt(SIGMA(2,2)));
fx1 = @(x) normpdf(x, MU(1), sqrt(SIGMA(1,1)));
fx2 = @(x) normpdf(x, MU(2), sqrt(SIGMA(2,2)));

% resistance
filename = ['fA_',num2str(t),'.mat'];

if ~exist(filename, 'file') || sim == 1
    construct_pdf(t);
end
load(filename, 'xA', 'fA')
load(['ffAA_',num2str(t),'.mat'],'xxAA','ffAA')
% fr = @(r) normpdf(r, 9000, 9000*0.1);
% fr = @(r) interp1(xA, fA, r);
fr = @(r) interp1(xxAA, ffAA, r);

% copula
switch copulatype
    
    case 'Gaussian'
        theta = ro_F;
%         theta = exp(-(delta_t/tau_F)^2);%sin(k_tau*pi/2);
%         fx1x2 = @(x1, x2) 1/sqrt(1-theta^2)*exp(-(norminv(Fx1(x1)).^2*theta^2 - 2*theta*norminv(Fx1(x1)).*norminv(Fx2(x2)) + norminv(Fx2(x2)).^2*theta^2)./(2*(1-theta^2))).*fx1(x1).*fx2(x2);
        fx1x2 = @(x1, x2) binormpdf(x1, x2, MU, SIGMA);
        %cpdf = @(x1, x2) reshape(copulapdf('Gaussian', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2); %for verification
    
    case 't' % dof = 2
        %theta = ro_F;
        theta = sin(k_tau*pi/2);
        %cpdf = @(x1, x2)
        fx1x2 = @(x1, x2) reshape(copulapdf('t', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta, 2), size(x1,1), []).*fx1(x1).*fx2(x2);
    
    case 'Clayton'
        theta = 2*k_tau/(1-k_tau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Clayton', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
   case 'rotClayton'
        theta = 2*k_tau/(1-k_tau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Clayton', [reshape(1-Fx1(x1),[],1), reshape(1-Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
        % not yet tested
        %fx1x2 = @(x1, x2) theta*(1+(theta+1)*(Fx1(x1)+Fx2(x2)-2*Fx1(x1)*Fx2(x2)))/((1+(theta-1)*(Fx1(x1)+Fx2(x2)))^2-4*Fx1(x1)*Fx2(x2)*theta*(theta-1))^(3/2);
    case 'Gumbel'
        theta = 1/(1-k_tau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Gumbel', [reshape(Fx1(x1),[],1), reshape(Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
    case 'rotGumbel'
        theta = 1/(1-k_tau);
        fx1x2 = @(x1, x2) reshape(copulapdf('Gumbel', [reshape(1-Fx1(x1),[],1), reshape(1-Fx2(x2),[],1)], theta), size(x1,1), []).*fx1(x1).*fx2(x2);
end

%==========================
% ANALYSIS/CALCULATION
%==========================
% integration limits
fx1max    = @(r) r;
fx2min    = @(r,x1) r;
frx1x2    = @(r,x1,x2) fr(r).*fx1x2(x1,x2);
tic
Pf_t = integral3(frx1x2, 500,9500 ,500,fx1max, fx2min,9500, 'AbsTol',1e-12, 'RelTol',1e-8);
tint = toc;
disp(['Numerical integration:   ', num2str(tint), ' seconds.']);
Pf_t
disp('=====================')

% Pf_d
nu_p(i) = Pf_t/delta_t;

% initial probability of failure, t=0
frx1 = @(r,x1) fr(r).*fx1(x1);
Pf_0 = integral2(frx1, 500,9500 ,fx1max,9500, 'AbsTol',1e-12, 'RelTol',1e-8)
end
disp('nu_p:')
nu_p'

end