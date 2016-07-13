clear all
close all
clc

% OPTIONS
%simulation number
N = 3;

L_meshed    = 1:1:1*365;
corr_length = 100;
pow         = 2;

%mean value
mu = 3500;

%standard deviation
sigma = 0.20*mu;


% CALCULATION
n_elem      = length(L_meshed)-1;
L_midpoint  = L_meshed(1:end-1) + diff(L_meshed)/2;


cov_matrix = element_cov_matrix(L_meshed, corr_length, pow);
%Cholesky decomposition of cov_matrix
% A = chol(cov_matrix,'lower');
% 
% ir1 = normrnd(mu, sigma, n_elem, N);
% dr1 = A*ir1;

r = mvnrnd(ones(n_elem,1)*mu, cov_matrix*sigma^2, N);

rr = zeros(size(r,1), length(L_meshed)*10);
calc_point = linspace(L_midpoint(1), L_midpoint(end), length(L_meshed)*10);
for i=1:size(r,1)
    rr(i,:) = interp1(L_midpoint, r(i,:), calc_point ,'spline');
end
%rr= r;
f_diag = figure('Position',[200, 200, 1200, 500]);
%f1 = plot(rr','black','LineWidth',2)   %black
f1 = plot(calc_point, rr','-','LineWidth',1);       %colored
grid on
hold on
% plot(calc_point, mean(rr),'blue','LineWidth',3)
% plot(calc_point, mean(rr)+std(rr),'red','LineWidth',3)
% plot(calc_point, mean(rr)-std(rr),'red','LineWidth',3)
ylim([0.9*min(min(rr)), 1.1*max(max(rr))])
% legend(f1,'simulation',...)
%         'Location','EastOutside')

% q = 500:100:7000;
% F = q;
% for i = 1:length(q)
%     id = max(rr,[],2) > q(i);
%     F(i) = sum(id)/N;
% end
% %save('p1.mat','p')
% qf          = q(1:end-1)+diff(q)/2;
% 
% f           = diff(F);
% scale       = trapz(qf,f);
% f           = f/scale;
% mu          = trapz(qf,f.*qf);
% var         = trapz(qf,f.*(qf-mu).^2);
% cov         = sqrt(var)/mu;
% 
% figure
% plot(q ,F)
% hold on
% plot(q, 1-normcdf(q, mu, mu*cov), 'r')
% grid on
% legend('numeric', 'mom normal')
% 
% mu
% cov

