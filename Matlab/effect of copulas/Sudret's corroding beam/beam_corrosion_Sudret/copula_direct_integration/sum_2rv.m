% Probability distribution function (pdf) of the sum of two independent random variables
%
% Convolution of two random variable -> performed using discrete Fourier transform
%
%SYNOPSYS:
% [fy, y] = SUM_2RV(pd1, pd2, n)
%
%INPUT:
% pd1
%  .fx_fun      handler of X1's pdf
%  .mean        mean value of X1
%  .min         'min' of X1 (to cover almost entirely the pdf's area
%  .max         'max' of X1 (to cover almost entirely the pdf's area 
%
% pd2
%OPTIONAL:
% n             number of point used in discrete Fourier transformation   
%
%OUTPUT:
% pdy
%
%NOTE(S):
% Y = X1 + X2
% the pdf of Y is shifted to the correct position based on the equivalence of mean values:
% E[Y] = E[X1 + X2]
% this requires that the rvs are given in a sufficently wide domain (covers almost all the densities)
% it would be great to be able to extract the shift without req of 'full' domain
%
% FAR FROM BEING GENERAL OR COMPLETE!

function pdy = sum_2rv(pd1, pd2, n, sign)
    
    % pre-processing
    if nargin < 4
        sign = 'add';
    end
    if nargin < 3
        n = 1e3;
    end
    
    switch lower(sign)
        case {'a','add','addition'}
            % do nothing
        case {'s','sub','subtract','subtraction'} %in case of sign = 'subtract' the order of rvs is important! X1-X2
            pd2.fx_fun  = @(x) pd2.fx_fun(-x);
            pd2.mean    = -pd2.mean;
            pd2.max     = -pd2.min;
            pd2.min     = -pd2.max;
        otherwise
            error(['Unknown ''sign'': ', sign])
    end
    
    % discretize the wider min-max range with n points - might not be sufficient for the smaller one!
    if abs(pd1.max - pd1.min) > abs(pd2.max - pd2.min)
        x1      = linspace(pd1.min, pd1.max, n);
        dx      = abs(diff(x1(1:2)));
        x2      = pd2.min:dx:pd2.max;
    else
        x2      = linspace(pd2.min, pd2.max, n);
        dx      = abs(diff(x2(1:2)));
        x1      = pd1.min:dx:pd1.max;        
    end

    % discrete points of the pdfs
    fx1     = pd1.fx_fun(x1);
    fx2     = pd2.fx_fun(x2);
    
    % convolution
    fy      = dx*conv(fx1, fx2, 'full');
    xx      = (0:numel(fy)-1)*dx;    

    % check the area under Y's pdf
    A       = trapz(xx,fy);
    if A < 1-1e-3
        warning(['The area under the pdf of the sum is considerably smaller than 1! A = ' num2str(A),...
            '. The [min, max] range or number of points (n) should be modified!'])
    end
    
    % shift to the correct position of Y's pdf
    mean_y  = trapz(xx,fy.*xx)/A;
    shift   = mean_y - (pd1.mean + pd2.mean);
    
    % function of Y's pdf
    fy_fun  = @(y) interp1(xx-shift, fy, y);
    
    % output /the min max might be refined, with small quantiles/
    pdy.mean    = mean_y-shift;
    pdy.min     = min(xx-shift);
    pdy.max     = max(xx-shift);
    pdy.fx_fun  = fy_fun;
    
end