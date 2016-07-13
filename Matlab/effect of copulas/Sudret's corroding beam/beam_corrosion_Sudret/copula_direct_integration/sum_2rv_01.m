% Probability distribution function (pdf) of the sum of two independent random variables
%
% Convolution of two random variable -> performed using discrete Fourier transform
%
%SYNOPSYS:
% [fy, y] = SUM_2RV(fx1, fx2)
%
%INPUT:
% fx1_fun       handler of X1's pdf
% fx2_fun       handler of X2's pdf
%
%
%OUTPUT:
%
%
%NOTE(S):
% Y = X1 + X2
% the pdf of Y is shifted to the correct position based on the equivalence of mean values:
% E[Y] = E[X1 + X2]
% this requires that the rvs are given in a sufficently wide domain (covers almost all the densities)
% it would be great to be able to extract the shift without req of 'full' domain
%
% FAR FROM BEING GENERAL OR COMPLETE!

function [fy_fun, y] = sum_2rv(fx1_fun, fx2_fun, n)
    
    % pre-processing
    if nargin < 3
        n = 1e3;
    end
    
%     % very inefficient, it would be better if the user supplied
%     Fx1_fun = @(x) integral(fx1_fun, -Inf, x); % WARNING!
%     x0      = fzero(@(x) Fx1_fun(x)-0.5, 1);
%     xmin    = fzero(@(x) Fx1_fun(x)-1e-8, x0);
%     xmax    = fzero(@(x) Fx1_fun(x)-(1-1e-8), x0);
    xmax = 8;
    xmin = -8;
    x1      = linspace(xmin, 2*xmax, n);
    x2      = x1; % WARNING!
    dx      = diff(x1(1:2));

    % discrete points of the pdf
    fx1     = fx1_fun(x1);
    fx2     = fx2_fun(x2);
    
    % mean of X1 and X2, needed for shifting Y's pdf
    mean_x1 = trapz(x1, fx1.*x1)/trapz(x1, fx1);
    mean_x2 = trapz(x2, fx2.*x2)/trapz(x2, fx2);
    

    % convolution
    fy      = dx*conv(fx1, fx2, 'same');
    y       = x1; % WARNING!

    % shift to the correct position of Y's pdf
    mean_y  = trapz(x1,fy.*x1)/trapz(x1,fy);
    shift   = mean_y - (mean_x1 + mean_x2);
    
    % function of Y's pdf
    fy_fun  = @(y) interp1(x1-shift, fy, y);
end