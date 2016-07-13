clear all
close all
clc

% OPTIONS
%simulation number
N = 1;

L_meshed    = 1:1:200;%*365;
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
A = chol(cov_matrix,'lower');

pon = 1;
a = mu - sqrt(3)*sigma;
b = mu + sqrt(3)*sigma;
sigma_ = sqrt(pon*sigma^2 + pon*((1-pon)*(a - mu))^2 + (1-pon)*pon*(mu-a));
mu_ = a*(1-pon) + pon*mu;
%ru = (-sqrt(3)*sigma) + rand(n_elem, N)*2*(sqrt(3)*sigma);
% ru = rand(n_elem, N);
% ion = ru <= pon;
% 
% ir1 = zeros(n_elem,1);
% ir1(~ion) = 0;
% ir1(ion) = (ru(ion) - (1-pon))*(b-a)/pon;
% 
% % 
% % ir1 = ru;%normrnd(0, sigma, n_elem, N);
% %ir1 = mvnrnd(ones(n_elem,1)*mu, eye(size(cov_matrix))*sigma^2, N)';%
% dr1 = A*ir1;
% r1 = dr1' + mu_;

r = mvnrnd(ones(n_elem,1)*mu_, cov_matrix*sigma^2, N);

rr = zeros(size(r,1), length(L_meshed)*10);
calc_point = linspace(L_midpoint(1), L_midpoint(end), length(L_meshed)*10);
for i=1:size(r,1)
    rr(i,:) = interp1(L_midpoint, r(i,:), calc_point ,'spline');
end
%rr= r;
f_diag = figure('Position',[200, 200, 1200, 500]);
%f1 = plot(rr','black','LineWidth',2)   %black
f1 = plot(calc_point, rr','-','LineWidth',1,'Color',[0.7,0.7,0.7]);       %colored
grid on
hold on
% plot(L_midpoint, r1, 'r')
% plot(calc_point, mean(rr),'blue','LineWidth',3)
% plot(calc_point, mean(rr)+std(rr),'red','LineWidth',3)
% plot(calc_point, mean(rr)-std(rr),'red','LineWidth',3)
%ylim([0.9*min(min(rr)), 1.1*max(max(rr))])
% legend(f1,'simulation',...)
%         'Location','EastOutside')

figure
autocorr(r,100)
% q = 3000:100:7000;
% p = q;
% for i = 1:length(q)
%     id = max(rr,[],2) > q(i);
%     p(i) = sum(id)/N;
% end
% %save('p1.mat','p')
% 
% figure
% plot(q,p)
% grid on