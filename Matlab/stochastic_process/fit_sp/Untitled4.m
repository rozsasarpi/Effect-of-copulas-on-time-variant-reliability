clearvars
clc
close all

% tic
% for ii = 1:1000
%     copulapdf('t',[0.5,0.4],0.4,3);
% end
% toc
% 
% tic
% copulapdf('t',repmat([0.5,0.4],1000,1),0.4,3);
% toc

% load('hr_tau_delta.mat', 'hr_tau', 'hr_delta')
% u = [0.8, 0.8];
% tau = 0.8;
% delta   = interp1(hr_tau, hr_delta, tau);
% bihr_copulapdf(u, delta)
% tic
% for ii = 1:1000
%     bitcopula([0.5,0.5],0.4,3);
% end
% toc

% 
% % rng(333)
% % n = 1;
% % ndim = 3;
% % COV = [1, 0.5, 0.7;
% %        0.5, 1, 0.5;
% %        0.7, 0.5, 1];
% %
% % x = rand(n,3);
% % u = normcdf(x,1,2);
% % copulapdf('Gaussian', u, COV)*prod(normpdf(x,repmat(1,1,ndim),repmat(2,1,ndim)))
% % % copulacdf('t', u, COV, 3)
% %
% % mvnpdf(x, repmat(1,1,ndim), COV*2^2)
% 
load('D:\Working folder\Matlab working folder\snow load\obs_database\obs_database_swe.mat')
load('D:\Working folder\Matlab working folder\snow load\obs_database\database_code_swe.mat')

point_ID = 2547; % Budapest
x = obs_database(2:end,point_ID);
i = find(x==9999, 1, 'last' )+1;
x = x(i:end);
x(x == -1) = NaN;

% % lim     = 20;
% % idx_p   = 1;
% % jj      = 1;
% % idx     = 1;
% % while idx1 < length(x)
% %     idx_p   = find(x(idx_p:end) > lim, 1);
% %     p       = x(idx_p);
% %     
% %     if jj > 1
% %         if idx_p > 20 + idx_peak(jj-1)
% %             peak(jj)        = p;
% %             idx_peak(jj)    = idx_p;
% %         end
% %     else
% %         
% %     end
% %     
% %     idx_p = idx_p + 1;
% % end
% 
% % mid_years = [1, round(182.625:365.25:rowNumber), rowNumber];
% %     max_a = zeros(length(mid_years)-1,1);
% % 
% %     for i = 1:length(mid_years)-1
% %         max_a(i) = max(obs_data(mid_years(i):mid_years(i+1)));
% %     end
% 
% %%
% i = 15;
% idx1 = (1+(i-1)*365);
% idx2 = (200+(i-1)*365);
% plot(x(idx1:idx2))