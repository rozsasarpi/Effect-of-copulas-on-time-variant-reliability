% Sudret et al. (2002). Comparison of methods for computing the probability of failure in time-variant reliability using the outcrossing approach
% Bruno Sudret (2008). Analytical derivation of the outcrossing rate in time-variant reliability problems. DOI:10.1080/15732470701270058
% simply supported beam subjected to stochastic load (F) and corrosion
%
% NOTE: the threedimensional integration is converted to a one-dimensinal one over a rectangle, this allows multiprecision integral!

clearvars
close all
clc

mp.Digits(100);
%==========================
% OPTIONS
%==========================
tt = mp('[0, 2, 5, 10, 15, 20]');
sim = 0;        % calculate the pdf of the R-Msw random variable?, 0 - no; 1 - yes
% tt = 0;

% copulatype = 'Gaussian';
% copulatype = 't';
% copulatype = 'rotClayton';
copulatype = 'Gumbel';
% copulatype = 'rotGumbel';


nu_p = zeros(numel(tt),1);
for i = 1:numel(tt)
    %time instant
    t       = tt(i);
    % correlation 'length' for F process
    tau_F   = mp('1/365');
    % span
    L       = mp('5');
    % time increment
    delta_t = tau_F*mp('0.001');
    % Pearsson correlation between 'adjactent' F realizations
    % Gaussian
    ro_F    = exp(-(delta_t/tau_F)^2);
    % Cauchy
%     ro_F    = (1+(delta_t/tau_F)^2).^(-2);
    % Kendall tau
    k_tau   = 2*asin(ro_F)/mp('pi');
    
    % stochastic load
    MU      = mp('[3500, 3500]')*L/4;
    SIGMA   = (L/4*mp('3500*0.2'))^2*ones(2,2);
    SIGMA(~eye(2)) = SIGMA(~eye(2))*ro_F;
   
    % integration limits
    int_l   = mp('500');
    int_u   = mp('9500');
%     int_l   = mp('300');
%     int_u   = mp('13500');
    
    %% direct integration solution
    % marginal cdfs and pdfs
    
    % resistance
    filename = ['ffBB_0.05_',num2str(double(t)),'.mat'];
    
    if ~exist(filename, 'file') || sim == 1
        construct_pdf(t);
    end
%     load(filename, 'xA', 'fA')
    load(['ffBB_0.05_',num2str(double(t)),'.mat'],'xxAA','ffAA')
    % fr = @(r) normpdf(r, 9000, 9000*0.1);
    % fr = @(r) interp1(xA, fA, r);
    fr = @(r) interp1(xxAA, ffAA, r);
    
    % copula
    switch copulatype
        case 'Gaussian'
            theta = ro_F;
            F     = @(r) binorm_rect_mp(int_l, r, r, int_u, MU, SIGMA).*fr(r);
        case 't' % dof = 2
            %theta = ro_F;
            theta = sin(k_tau*pi/2);
            F     = @(r) bit_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta, 2).*fr(r);
        case 'Clayton'
            theta = 2*k_tau/(1-k_tau);
            error(['Not implemented copula:', copulatype])
        case 'rotClayton'
            theta = 2*k_tau/(1-k_tau);
            F     = @(r) birotclay_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta).*fr(r);
        case 'Gumbel'
            theta = 1/(1-k_tau);
            F     = @(r) bigumb_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta).*fr(r);
        case 'rotGumbel'
            theta = 1/(1-k_tau);
            F     = @(r) birotgumb_rect_mp(int_l, r, r, int_u, MU, sqrt([SIGMA(1), SIGMA(4)]), theta).*fr(r);
        otherwise
            error(['Unknown copula type:', copulatype])
    end
    
    %==========================
    % ANALYSIS/CALCULATION
    %==========================
    tic;
    Pf_t = quadgk(F, int_l, int_u, 'RelTol',mp('1e-8'),'AbsTol',mp('1e-12'))
%     Pf_t = quadgk(F, int_l, int_u)
    toc;
    tint = toc;
    disp(['Numerical integration:   ', num2str(tint), ' seconds.']);
    disp('=====================')
    
    % Pf_d
    nu_p(i) = Pf_t/delta_t;
    
end
disp('nu_p:')
nu_p'