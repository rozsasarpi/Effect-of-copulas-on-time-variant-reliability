% clear all
% close all
% clc

%CUSTOM FUNCTION(S):
%sum_2rv.m
%prod_2rv.m

%[mm], [year], [N]

function f = construct_pdf(t, x, method, fig)
    
    if nargin < 4
        fig = 'no';
    end
    
    %=======================
    % OPTIONS
    %=======================
    d_scale = 1e-3; %if d_scale = 1 -> mm
    qq      = 1e-8;
    
    %=======================
    % PARAMETERS
    %=======================
    % elapsed time
    % t       = 10;
    % span of the beam
    L       = 5000*d_scale;
    % corrosion rate
    kappa   = 0.05*d_scale;   % mm/year
    % unit weight of the beam [N/mm^3]
    ro      = 78.5*10^-6*d_scale^-3;
    % reduction of dimensions
    dd    = 2*kappa*t;
    
    % yield strength
    fy_mean         = 240*d_scale^-2;
    fy_cov          = 0.10;
    [fy_mu, fy_sigma] = lognormpar(fy_mean, fy_cov);
    
    % beam width
    b0_mean         = 200*d_scale;
    b0_cov          = 0.05;
    [b0_mu, b0_sigma] = lognormpar(b0_mean, b0_cov);
    
    % beam heigth
    h0_mean         = 40*d_scale;
    h0_cov          = 0.10;
    [h0_mu, h0_sigma] = lognormpar(h0_mean, h0_cov);
    
    switch lower(method)
        case 'simulation'
            %% Simulation
            N = 1e6;
            n = 100;
            tic
            sampleA = zeros(N,n);
            parfor i = 1:n
                fy  = lognrnd(fy_mu, fy_sigma, N, 1);
                b0  = lognrnd(b0_mu, b0_sigma, N, 1);
                h0  = lognrnd(h0_mu, h0_sigma, N, 1);
                
                b       = b0 - dd;
                h       = h0 - dd;
                sA      = b.*h.^2.*fy./4 - ro.*b0.*h0.*L.^2/8;
                
                sampleA(:,i) = sA;
            end
            sampleA = sampleA(:);
            tsim = toc;
            disp(['Simulation:              ', num2str(tsim), ' seconds.']);
            % load(['tmp\sample_',num2str(ID),'.mat'], 'sample')
            %kernel density implementation from file exchange (Botev)
            % with default 2^14 meshpoints and min  max (range) +0.1*range
            tic
            %     sampleA = sampleA(sampleA>2000 & sampleA<9000);
            %     numel(sampleA)
            %     numel(sampleA)/(N*n)
            [~,fA,xA] = kde(sampleA);
            tkde = toc;
            disp(['Kernel construction:     ', num2str(tkde), ' seconds.']);
            
            save(['fB_',num2str(t),'.mat'], 'xA', 'fA')
            
            %PLOT
            switch lower(fig)
                case 'yes'
                    [f, x] = hist(sampleA,1000);
                    f = f/trapz(x,f);
                    
                    plot(x,f)
                    hold on
                    plot(xA,fA,'r')
                    grid on
                    legend('norm hist','kernel')
                case 'no'
                    %do nothing
            end
            cfd
            %% numerical integration
        case 'integration'
            x   = x*d_scale;
            h   = (eps)^(1/3)*max(1, abs(x));
            xph = x + h;
            xmh = x - h;
            dx  = xph - xmh; % to account for rounding error
            
            a   = [xph, xmh];
            Fn  = zeros(2,1);
%             tic
            for i = 1:2
                
                ffy    = @(x) lognormpdf(x, fy_mean, fy_cov);
                fb0    = @(x) lognormpdf(x, b0_mean, b0_cov);
                fh0    = @(x) lognormpdf(x, h0_mean, h0_cov);
                % integration limits
                ffymax    = @(b0,h0) 4./((b0-dd).*(h0-dd).^2).*(a(i) + ro.*b0.*h0.*L.^2/8);
                fb0h0fy   = @(b0,h0,fy) fb0(b0).*fh0(h0).*ffy(fy);
                
                Pf_d = integral3(fb0h0fy, 100*d_scale,300*d_scale ,20*d_scale,60*d_scale, 50*d_scale^-2,ffymax, 'AbsTol',1e-33, 'RelTol',1e-5);
                Fn(i) = Pf_d;
%                 disp(Pf_d)
            end
%             toc
            
            f = (Fn(1) - Fn(2))./(dx);
%             disp(f)
        case 'rv algebra'
            
            %% Sum and product of distributions - would be more elegent if it would work..
            % %=======================
            % % PREPARATION of stucture variables for convolution
            % %=======================
            % % fy
            % pd_fy.mean      = fy_mean;
            % pd_fy.min       = 0;
            % pd_fy.max       = fzero(@(x) lognormcdf(x, fy_mean, fy_cov) - (1-qq), fy_mean);
            % pd_fy.fx_fun    = @(x) lognormpdf(x, fy_mean, fy_cov);
            %
            % % b
            % pd_b0.mean      = b0_mean;
            % pd_b0.min       = 0;
            % pd_b0.max       = fzero(@(x) lognormcdf(x, b0_mean, b0_cov) - (1-qq), b0_mean);
            % pd_b0.fx_fun    = @(x) lognormpdf(x, b0_mean, b0_cov);
            %
            % % h
            % pd_h0.mean      = h0_mean;
            % pd_h0.min       = 0;
            % pd_h0.max       = fzero(@(x) lognormcdf(x, h0_mean, h0_cov) - (1-qq), h0_mean);
            % pd_h0.fx_fun    = @(x) lognormpdf(x, h0_mean, h0_cov);
            %
            % tic
            % %=======================
            % % CALCULATION - constructing the pdf
            % %=======================
            % % PDF of the moment from the self-weigth of the beam
            % sigma           = sqrt(b0_sigma^2 + h0_sigma^2);
            % mu              = b0_mu + h0_mu;
            % [M, V]          = lognstat(mu, sigma);
            % pd_msw.mean     = M;
            % pd_msw.min      = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - qq, M);
            % pd_msw.max      = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - (1-qq), M);
            % pd_msw.fx_fun   = @(x) lognormpdf(x, M, sqrt(V)/M);
            %
            % pd_Msw.mean     = ro*L^2/8*pd_msw.mean;
            % pd_Msw.min      = ro*L^2/8*pd_msw.min;
            % pd_Msw.max      = ro*L^2/8*pd_msw.max;
            % pd_Msw.fx_fun   = @(x) lognormpdf(x, ro*L^2/8*M, sqrt(ro*L^2/8*V)/(ro*L^2/8*M));
            %
            % % PDF of the plastic moment resistance of the beam
            % % 1/4*[b0*h0^2 - 2*dd*b0*h0 + dd^2*b0 - dd*h0^2 + 2*dd^2*h0 (- dd^3)]*fy
            %
            % % basic terms 1x
            % %11 _ b0*h0^2
            % sigma           = sqrt((b0_sigma)^2 + 2*(h0_sigma)^2);
            % mu              = b0_mu + 2*h0_mu;
            % [M, V]          = lognstat(mu, sigma);
            % pd_11.mean      = M;
            % pd_11.min       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - qq, M);
            % pd_11.max       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - (1-qq), M);
            % pd_11.fx_fun    = @(x) lognormpdf(x, M, sqrt(V)/M);
            %
            % %12 _ 2*dd*b0*h0 /minus sign added later/
            % sigma           = sqrt((b0_sigma)^2 + (h0_sigma)^2);
            % mu              = b0_mu + h0_mu;
            % [M, V]          = lognstat(mu, sigma);
            % M               = 2*dd*M;
            % V               = (2*dd)^2*V;
            % pd_12.mean      = M;
            % pd_12.min       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - qq, M);
            % pd_12.max       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - (1-qq), M);
            % pd_12.fx_fun    = @(x) lognormpdf(x, M, sqrt(V)/M);
            %
            % %13 _ dd^2*b0
            % sigma           = b0_sigma;
            % mu              = b0_mu;
            % [M, V]          = lognstat(mu, sigma);
            % M               = dd^2*M;
            % V               = (dd^2)^2*V;
            % pd_13.mean      = M;
            % pd_13.min       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - qq, M);
            % pd_13.max       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - (1-qq), M);
            % pd_13.fx_fun    = @(x) lognormpdf(x, M, sqrt(V)/M);
            %
            % %14 _ dd*h0^2 /minus sign added later/
            % sigma           = sqrt(2)*h0_sigma;
            % mu              = 2*h0_mu;
            % [M, V]          = lognstat(mu, sigma);
            % M               = dd*M;
            % V               = dd^2*V;
            % pd_14.mean      = M;
            % pd_14.min       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - qq, M);
            % pd_14.max       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - (1-qq), M);
            % pd_14.fx_fun    = @(x) lognormpdf(x, M, sqrt(V)/M);
            %
            % %15 _ 2*dd^2*h0
            % sigma           = h0_sigma;
            % mu              = h0_mu;
            % [M, V]          = lognstat(mu, sigma);
            % M               = 2*dd^2*M;
            % V               = (2*dd^2)^2*V;
            % pd_15.mean      = M;
            % pd_15.min       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - qq, M);
            % pd_15.max       = fzero(@(x) lognormcdf(x, M, sqrt(V)/M) - (1-qq), M);
            % pd_15.fx_fun    = @(x) lognormpdf(x, M, sqrt(V)/M);
            %
            % % convolution 2x
            % pd_21 = sum_2rv(pd_11, pd_12, 1e5, 'subtraction');
            % pd_22 = sum_2rv(pd_13, pd_15, 1e5);
            % pd_23 = sum_2rv(pd_21, pd_14, 1e6, 'subtraction');
            % pd_24 = sum_2rv(pd_23, pd_22, 1e6);
            %
            % % shifting with -dd^3
            % pd_Wpl.mean      = pd_24.mean - dd^3;
            % pd_Wpl.max       = pd_24.max - dd^3;
            % pd_Wpl.min       = pd_24.min - dd^3;
            % pd_Wpl.fx_fun    = @(x) pd_24.fx_fun(x+dd^3);
            %
            % % product distribution
            % % pd31.fx_fun*fy = Wpl*fy
            % pd_r = prod_2rv(pd_Wpl, pd_fy);
            % pd_r.mean        = pd_Wpl.mean*pd_fy.mean;
            % pd_r.max         = fzero(@(x) pd_r.fx_fun - qq, M);
            %
            % % R = 1/4*Wpl*fy
            % pd_R.mean      = 1/4*pd_r.mean;
            % pd_R.max       = 1/4*pd_r.max;
            % pd_R.min       = 1/4*pd_r.min;
            % pd_R.fR_fun = @(x) pd_r.fx_fun(4*x);
            % toc
            %
            % % final distribution to be used in further analysis
            % % A = R - Msw
            % pd_A = sum_2rv(pd_R, pd_Msw, 1e6, 'substraction');
            
            %R and Msw are correlated... and h0^2 should be considered as correlated..
        otherwise
            error(['Unknown method:', method])
    end
    
    function [mu_lognorm, sigma_lognorm] = lognormpar(meanX, covX)
        mu_lognorm      = log(meanX.^2./sqrt((meanX.*covX).^2 + meanX.^2));
        sigma_lognorm   = sqrt(log((meanX.*covX).^2./meanX.^2 + 1));
    end
end



