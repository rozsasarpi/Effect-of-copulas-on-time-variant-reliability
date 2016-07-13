% Shanks transformation to increase rate of convergence of limit

function shanks_trans_test
% clear variables
close all
clc

% https://en.wikipedia.org/wiki/Shanks_transformation#Example
% k and n are exchanged
% fun     = @(n) 4*sum((-1).^(0:n).*(1./(2*(0:n)+1)));
ffun    = @(x) sin(x)./x;
fun     = @(k) ffun(2.^(-k));
% fun     = @(k) ffun(1./(k));
% fun     = @(k) ffun(1./factorial(k));

lim_x   = 1;

k       = 5;

ini  = 0;
format longE
A    = arrayfun(fun, ini:(ini+k+1)).';
SA1  = shanks_trans(A);
SA2  = shanks_trans(SA1);
SA3  = shanks_trans(SA2);

[A, SA1, SA2, SA3]

[A, SA1, SA2, SA3] - lim_x

min(2.^(-(ini:(ini+k+1))))
% min(1./(0:(k+1)))
% min(1./factorial(k))

% [A(1:(end-1))-lim_x, SA-lim_x]

% apply Richardson extrapolation to the Shanks transformed sequence
% kk      = 1:k;
% seq     = (2.^(-kk)).';
% p       = polyfit(seq, SA(2:end), k-1);
% 
% pf      = polyval(p, seq)-1
% 
% polyval(p, 0)-1

    function SA = shanks_trans(A)        
        nn   = length(A);
        SA   = nan(nn,1);
        for ii = 0:(nn-1)
            idx = ii + 1;
            if ii == 0
                SA(idx) = NaN;
            elseif ii > nn-2
                SA(idx) = NaN; 
            else
                SA(idx) = (A(idx+1).*A(idx-1) - A(idx).^2)./(A(idx+1) - 2*A(idx) + A(idx-1));
            end
        end
    end


end