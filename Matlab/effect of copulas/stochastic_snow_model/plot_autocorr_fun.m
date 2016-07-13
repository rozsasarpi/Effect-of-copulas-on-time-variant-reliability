% plot autocorrelation function

clearvars
close all
clc

x = 0:0.1:6;
xx = sort([-x, x]);

Corr.length = 2;
Corr.pow    = 2;

figure('Position', [400,100,700,280])
% AUTOCORR FUN PLOT
% subplot(1,2,1)
% exponential
Corr.type   = 'exp';
CM = cov_matrix(x, Corr);
% plot(xx, [fliplr(CM(1,:)), CM(1,:)], 'black', 'LineStyle', '-.', 'LineWidth', 1)

% % chauchy
hold on
Corr.type   = 'c';
CM = cov_matrix(x, Corr);
plot(xx, [fliplr(CM(1,:)), CM(1,:)], 'black', 'LineStyle', '-', 'LineWidth', 1)

% % gaussian
Corr.type   = 'g';
CM = cov_matrix(x, Corr);
plot(xx, [fliplr(CM(1,:)), CM(1,:)], 'black', 'LineStyle', '--', 'LineWidth', 1)

xmin = min(xlim);
% plot([xmin, xmin], [0, 1], 'white', 'LineWidth', 2)

% legend('exponential', 'cauchy', 'gaussian')
xlabel('$\Delta t$', 'Interpreter', 'LaTeX')
ylabel('$\rho$', 'Interpreter', 'LaTeX')
box on
hl = legend('Cauchy', 'Gauss', 'Location', 'EastOutside');
% hl.Interprete = 'LaTeX';
hl.Box = 'off';
set(gca,'XTick', [-6:2:6])

% % EQUATIONS
% subplot(1,2,2)
% % plot([0, 0.15], [0.8, 0.8], 'black', 'LineStyle', '-.', 'LineWidth', 1)
% hold on
% plot([0, 0.15], [0.6, 0.6], 'black', 'LineStyle', '-', 'LineWidth', 1)
% hold on
% plot([0, 0.15], [0.3, 0.3], 'black', 'LineStyle', '--', 'LineWidth', 1)
% xlim([0,1])
% ylim([0,1])
% % text(0.2, 0.8, 'Exponential:  $\rho (\Delta t) = \exp \left( { - \frac{{\Delta t}}{{{\tau _F}}}} \right)$', 'Interpreter', 'LaTeX')
% text(0.2, 0.6, 'Cauchy:  $\rho (\Delta t) = {\left( {1 + {{\left( {\frac{{\Delta t}}{{{\tau _F}}}} \right)}^2}} \right)^{ - 2}}$', 'Interpreter', 'LaTeX')
% text(0.2, 0.3, 'Gaussian:  $\rho (\Delta t) = \exp \left( { - {{\left( {\frac{{\Delta t}}{{{\tau _F}}}} \right)}^2}} \right)$', 'Interpreter', 'LaTeX')
% axis off
% box off