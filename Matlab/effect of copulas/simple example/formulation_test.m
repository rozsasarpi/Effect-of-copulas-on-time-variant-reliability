clearvars
close all
clc

R   = mp('3');
rho = mp('0.999');
% integral2(@(x,y) mvnpdf([x, y], [0,0], [1, 0.5; 0.5, 1]), -Inf,R, R,Inf)
integral2(@(x,y) reshape(binormpdf([x(:), y(:)], [0,0], double([1, rho; rho, 1])), size(x)), -Inf,double(R), double(R),Inf, 'Abstol', 1e-16, 'Reltol', 1e-16)

p_s = normcdf(double(R));

% prone to catastrophic cancelation!!
p_s - copulacdf('Gaussian', [p_s, p_s], double(rho))

u = repmat(1-bivnor_mp(R,mp('-Inf'),0,50),1,2);

p_s - binorm_copulacdf_mp(u, rho)
