%Generates covariance matrix for 1D Gaussian distributed (spatially) random field
%
% The random field is discretized in the midpoint of pfem finite elements.
%
%SYNOPSIS
% cov_matrix = ELEMENT_COV_MATRIX(L_meshed, stdev, corr_length, pow)
%
%INPUT
% L_meshed      - coordinate of the pfem mesh nodes in the element local
%                 coordinate system /vector, nx1/
% corr_length   - correlation length of the random variable/field /constant/
%
%OPTIONAL
% pow           - power/exponent (larger means faster decrease in
%                 correlation with distance), /constant (2)/
%
%OUTPUT
% cov_matrix    - covariance matrix of the random field in the mid points
%                 of pfem finite elements
%                 /matrix, nxn/
%COMMENTS
% (1) NOTE: The script uses Statistical Toolbox's functions (pdist, squareform).

function cov_matrix = element_cov_matrix(L_meshed, corr_length, pow)

    if isrow(L_meshed) == 1
        L_meshed = L_meshed';
    end

    if nargin < 4
        pow = 2;    % default value for the power
    end
    
    % midpoints of the pfem finite elements
    L_midpoint = L_meshed(1:end-1) + diff(L_meshed)/2;

    % Eucledian distance of the pfem mesh nodes
    dx = squareform(pdist(L_midpoint,'euclidean'));

    % covariance matrix
    cov_matrix = 1^2.*exp(-dx./corr_length).^pow;

end