clearvars
close all
clc
format long


%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% X1 and X2 mean and variance-covariance matrix
MU      = [3,2];
SIGMA   = [1, 0.5; 0.5, 2];
mu1     = MU(1);
mu2     = MU(2);
sx1     = sqrt(SIGMA(1));
sx2     = sqrt(SIGMA(4));
rho     = SIGMA(2)/sqrt(SIGMA(1)*SIGMA(4));

% limits for testing
XL      = [1,2];
XU      = [2,3];


%--------------------------------------------------------------------------
% DISTRIBUTIONS
%--------------------------------------------------------------------------
fx1x2   = @(x1,x2) reshape(mvnpdf([x1(:), x2(:)], MU, SIGMA), size(x1));

% Fx1x2   = @(xl, xu) mvncdf(xu, mu, SIGMA) - mvncdf([xu(1),xl(2)], mu, SIGMA) -...
%     mvncdf([xl(1),xu(2)], mu, SIGMA) + mvncdf(xl, mu, SIGMA);
% mp fucks up the calculations, thus double is needed...
% Fx1x2   = @(xl, xu) double(bivnor(-(xu(1)-mu(1))/sx1, -(xu(2)-mu(2))/sx2, rho) - bivnor(-(xu(1)-mu(1))/sx1, -(xl(2)-mu(2))/sx2, rho) -...
%     bivnor(-(xl(1)-mu(1))/sx1, -(xu(2)-mu(2))/sx2, rho) + bivnor(-(xl(1)-mu(1))/sx1, -(xl(2)-mu(2))/sx2, rho));

Fx1x2   = @(xl, xu) bivnor(-(xu(1)-MU(1))/sx1, -(xu(2)-MU(2))/sx2, rho) - bivnor(-(xu(1)-MU(1))/sx1, -(xl(2)-MU(2))/sx2, rho) -...
    bivnor(-(xl(1)-MU(1))/sx1, -(xu(2)-MU(2))/sx2, rho) + bivnor(-(xl(1)-MU(1))/sx1, -(xl(2)-MU(2))/sx2, rho);

% binormpdf_rect_mp(XL(1), XL(2), XU(1), XU(2), MU, SIGMA)
Fx1x2(XL,XU)
mvncdf(XL, XU, MU, SIGMA)

fr      = @(r) normpdf(r, 0, 1);

% 3D EXAMPLE
%--------------------------------------------------------------------------
% SIMULATION
%--------------------------------------------------------------------------
N       = 1e6;
rx1x2   = mvnrnd(MU,SIGMA,N);
rr      = normrnd(0,1,N,1);
nf      = sum(rr-rx1x2(:,1) >= 0 & rr-rx1x2(:,2) < 0);
Pf_3D_MC = nf/N
%--------------------------------------------------------------------------
% INTEGRAL
%--------------------------------------------------------------------------

disp('========== 3D integration ==========')
fx1max  = @(r) r;
fx2min  = @(r,x1) r;
frx1x2  = @(r,x1,x2) fr(r).*fx1x2(x1,x2);

int_l   = -4;
int_u   = 6;
tic
Pf_3D_3 = integral3(frx1x2, int_l,int_u ,int_l,fx1max, fx2min,int_u, 'AbsTol',1e-10, 'RelTol',1e-06)
toc

disp('========== 1D integration ==========')
% F       = @(r) Fx1x2([int_l, r], [r, int_u]).*fr(r);
F       = @(r) binorm_rect_mp(int_l, r, r, int_u, MU, SIGMA).*fr(r);

tic
% Pf_3D_1 = integral(F, int_l,int_u, 'ArrayValued', 1, 'AbsTol',1e-10, 'RelTol',1e-06)
% Pf_3D_1 = quad(F, int_l,int_u, mp('1e-08'))
Pf_3D_1 = quadgk(F, int_l,int_u, 'RelTol',mp('1e-6'),'AbsTol',mp('1e-10'))
toc
% although it is quite time consuming this formulation allows multiprecision calculation!

% 2D EXAMPLE
% fx1 = @(x1) normpdf(x1, mu(1), SIGMA(1));
% Fx1 = @(il) normcdf(il(2), mu(1), SIGMA(1)) - normcdf(il(1), mu(1), SIGMA(1));
% 
% disp('========== 2D integration ==========')
% fx1max  = @(r) r;
% frx1    = @(r,x1) fr(r).*fx1(x1);
% 
% int_l   = -4;
% int_u   = 6;
% tic
% Pf_2D_2 = integral2(frx1, int_l,int_u ,int_l,fx1max, 'AbsTol',1e-10, 'RelTol',1e-06)
% toc
% 
% disp('========== 1D integration ==========')
% F       = @(r) Fx1([int_l, r]).*fr(r);
% 
% tic
% Pf_2D_1 = integral(F, int_l, int_u, 'ArrayValued', 1, 'AbsTol',1e-10, 'RelTol',1e-06)
% toc



