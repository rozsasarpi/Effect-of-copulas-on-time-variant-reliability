% Convolution RG = R-G
%
% [Probvar, pdRG] = CONV_RG(Probvar)

function [Probvar, pdRG] = conv_RG(Probvar)

rv  = {'R', 'G'};
nrv = length(rv);

pds(nrv).min = NaN;
% loop over random variables
for ii = 1:nrv
    
    if Probvar.(rv{ii}).dist == 1
        pd.fx_fun = @(x) normpdf(x, Probvar.(rv{ii}).mean, Probvar.(rv{ii}).cov);
        pd.min = Probvar.(rv{ii}).mean*(1-30*Probvar.(rv{ii}).cov);
    elseif Probvar.R.dist == 2
        pd.fx_fun = @(x) lognormpdf(x, Probvar.(rv{ii}).mean, Probvar.(rv{ii}).cov);
        pd.min = max(Probvar.(rv{ii}).mean*(1-30*Probvar.(rv{ii}).cov),0);
    else
        error('Not yet implemented distribution type')
    end
    
    pd.mean = Probvar.(rv{ii}).mean;
    pd.max  = Probvar.(rv{ii}).mean*(1+30*Probvar.(rv{ii}).cov);
    
    if ii == 1
        pdR = pd;
    elseif ii == 2
        pdG = pd;
    end
end

% CONVOLUTION, R-G
pdRG = sum_2rv(pdR, pdG, 1e5, 's');

% Probvar.RG should be added

end
