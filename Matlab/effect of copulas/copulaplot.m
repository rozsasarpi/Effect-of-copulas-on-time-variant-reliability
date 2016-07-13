% plots of copulas

clear all
close all
clc

%OPTIONS
%================================================
% Pearson correlation coefficient
rho = [0.1, 0.5, 0.9];
% ??
alpha = [0.7, 2, 3];
copulatype = 'Gaussian';
%================================================

switch lower(copulatype)
    case 'gaussian'
        param = rho;
        paramname = '\rho';
    case {'clayton', 'frank', 'gumbel'}
        param = alpha;
        paramname = '\alpha';
end

% PLOT
hFig = figure;
set(hFig, 'Position', [0 400 1200 500])
for i = 1:length(param)
    %grid points
    x = linspace(1e-2, 1-1e-2, 20);
    y = x;
    [X, Y] = meshgrid(x,y);
    % uniform marginals
    U = [X(:), Y(:)];
    
    % Gaussian copula
    Z = copulapdf(copulatype, U, param(i));
    Z = reshape(Z, length(x), []);
    
    subplot(1,length(rho),i)
    surf(X,Y,Z, 'FaceAlpha', 0.6)
    %zlim([min(min(Z)),min(max(max(Z)),4)])
    zlim([min(min(Z)),4])
    %colormap(summer)
    caxis([min(min(Z)),4])
    %     bg_color = get(gca,'Color');
    %     set(gca,'ZColor',bg_color,'ZTick',[])
    view([-20,40])
    xlabel('U_1')
    ylabel('U_2')
    zlabel('Probability Density')
    title([copulatype, ' copula, ', paramname,' = ',num2str(param(i))])
end

