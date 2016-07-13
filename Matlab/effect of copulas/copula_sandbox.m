
x = linspace(-3, 3, 20);
y = x;
[X, Y] = meshgrid(x,y);
% uniform marginals
F = [normcdf(X(:),0,1), normcdf(Y(:),0,1)];

Z = copulacdf('Gaussian', F, 0.1);
Z = reshape(Z, length(x), []);

surf(X,Y,Z, 'FaceAlpha', 0.6, 'FaceColor', 'blue')

NZ = mvncdf([X(:), Y(:)], [0, 0], [1, 0.1; 0.1, 1]);
NZ = reshape(NZ, length(x), []);

hold on
surf(X,Y,Z, 'FaceAlpha', 0.6, 'FaceColor', 'red')

%%
% G2    bivariate cdf of the copula
% g2    bivariate pdf of the copula
% G1    univariate cdf of the copula
% invF  univariate marginal cdf of the examined distribution
% all functions must accept arrays and output array, element-wise processing is required
clear all
close all
clc

mean_xi = 1;
mean_xj = 1;
stdv_xi = 1;
stdv_xj = 1;

invFi = @(P) lognorminv(P, mean_xi, stdv_xi/mean_xi);
invFj = @(P) lognorminv(P, mean_xj, stdv_xj/mean_xj);

ymin = -8;
ymax = 8;
rho_ij = 0.1;

copulatype = 'Gaussian';

switch lower(copulatype)
    case 'gaussian'
        % standard normal with correlation rho_e
        %g2 = @(yi, yj, rho_e) reshape(mvnpdf([yi(:), yj(:)], [0, 0], [1, rho_e; rho_e, 1]), size(yi,1), size(yi,2));
        g2 = @(yi, yj, rho_e) 1/(2*pi*sqrt(1-rho_e^2)) * exp(-1/(2*(1-rho_e^2)) * (yi.^2 - 2*rho_e*yi.*yj + yj.^2));
        G1 = @(y) normcdf(y);        
end

I = @(yi, yj, rho_e) (invFi(G1(yi)) - mean_xi).*(invFj(G1(yj)) - mean_xj).*g2(yi, yj, rho_e);

rho = @(rho_e) 1/(stdv_xi*stdv_xj)*integral2(@(yi, yj) I(yi, yj, rho_e), ymin, ymax, ymin, ymax);

rho_diff = @(rho_e) abs(rho(rho_e) - rho_ij);

tic
[x,fval,exitflag,output] = fmincon(rho_diff, rho_ij, [1;-1], [1;0]);%fzero(rho_diff, [0,1-1e-5]);
toc
disp(['The equivalent correlation coefficient: ', num2str(x)])