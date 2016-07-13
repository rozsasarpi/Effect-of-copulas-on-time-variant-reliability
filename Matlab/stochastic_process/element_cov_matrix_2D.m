%Generates covariance matrix for 2D Gaussian distributed (spatially) random field
%
% The random field is discretized in the midpoint of pfem finite elements.
%
%SYNOPSIS
% cov_matrix = ELEMENT_COV_MATRIX_2D(L_meshed, stdev, corr_length, pow)
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

function cov_matrix = element_cov_matrix_2D(Coord, corr_length, pow)

    if nargin < 4
        pow = 2;    % default value for the power
    end

    % Eucledian distance of the nodes
    dr = squareform(pdist(Coord,'euclidean'));

    % covariance matrix
    cov_matrix = 1^2.*exp(-dr./corr_length).^pow;

end