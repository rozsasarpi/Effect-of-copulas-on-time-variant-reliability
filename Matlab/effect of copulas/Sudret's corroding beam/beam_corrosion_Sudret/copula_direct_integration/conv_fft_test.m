clear all
close all
clc

% Y = X1 + X2
% the pdf of Y is shifted to the corrected position based on the equivalence of means
% E[Y] = E[X1 + X2]
% this requires that the rvs are given in a sufficently wide domain (covers the almost all the density) 
% it would be great to be able to extract the shift without the req of 'full' domain 

x1  = linspace(-10,10,1e4);
dx  = (x1(2)-x1(1));
x2  = -8:dx:8;

y1  = normpdf(x1,1,1);
y2  = normpdf(x2,2,0.1);
% y2  = normpdf(x2,1,1);
% y   = lognormpdf(x,1,1);

% n = 1e5;
% yr1 = lognormrnd(1,1,n,1);
% yr2 = lognormrnd(1,1,n,1);
% [ff, xx] = hist(yr1 + yr2, length(x));
% y  = rectangularPulse(1,2,x);
tic
mean_y1 = trapz(x1,y1.*x1)/trapz(x1,y1);
mean_y2 = trapz(x2,y2.*x2)/trapz(x2,y2);


% cy = dx*conv(y1,y2, 'same');
toc
cy      = dx*conv(y1,y2, 'full');
cx      = (0:numel(cy)-1)*dx;
mean_cy = trapz(cx,cy.*cx)/trapz(cx,cy);
shift   = mean_cy - (mean_y1 + mean_y2);

%----------
tic
pd1.mean = 1;
pd1.max = 10;
pd1.min = -10;
pd1.fx_fun = @(x) normpdf(x,1,1);

pd2.mean = 2;
pd2.max = 8;
pd2.min = -8;
pd2.fx_fun = @(x) normpdf(x,2,0.1);

pdy = sum_2rv(pd1, pd2, 1e3, 'a');
toc
%----------

ye = normpdf(cx-shift,3,sqrt(1+0.1^2));
% ye = normpdf(cx-shift,2,sqrt(2));
% ye = ff/trapz(xx,ff);

plot(x1,y1,'m')
hold on
plot(x2,y2,'m--')

plot(cx-shift,cy,'r-')
plot(cx-shift,ye,'--')
%----------
plot(x1,pdy.fx_fun(x1),'y--')
%----------
legend('Y1', 'Y2', 'Y+Y FFT', 'Y+Y exact')
ylim([0,1])

figure
dy = cy-ye;
plot(cx-shift, dy./ye*100, 'LineWidth', 2)
title('relative difference of FFT and exact solution [%]')
xlabel('x')
ylabel('(FFT-exact)/exact [%]')
grid on
line([min(cx-shift),max(cx-shift)], [0,0], 'Color', 'red')

disp(['abs(max(rel error) : ', num2str(max(abs(dy)./ye*100)), '%'])
