% Select annual maxima from daily observations
%
% [max_a, f_scatter] = ANNUAL_MAX(obs_data, rowNumber, years, plotFig)
%
% data      - daily observations
% rowNumber - number of observations
% years     - year of maxima (vector)
% plotFig   - plot? y(1),n(2); n(2) - default (str,integer)

function [max_a, max_idx, f_scatter] = annual_max(obs_data, rowNumber, years, plotFig)

    if nargin < 4
        plotFig = 'n';
        f_scatter = 0;
    else
        if isnumeric(plotFig) == 1
            if plotFig == 1
                plotFig = 'y';
            else
                plotFig = 'n';
            end
        end
    end

    mid_years   = [1, round(182.625:365.25:rowNumber), rowNumber];
    max_a       = zeros(length(mid_years)-1,1);
    max_idx     = zeros(length(mid_years)-1,1);

    for i = 1:length(mid_years)-1
        [max_a(i), idx] = max(obs_data(mid_years(i):mid_years(i+1)));
        max_idx(i) = mid_years(i) + (idx-1);
    end
    
    switch lower(plotFig)
        case 'y'
            f_scatter = figure('Position',[400, 400, 800, 400]);
            plot(years,max_a,'o-',...
                              'Color','black',...
                              'MarkerFaceColor', [0.7, 0.7, 0.7],...
                              'MarkerSize',5);
            p = polyfit(years,max_a,1);
            linFit = polyval(p,years);
            hold on
            plot(years,linFit, 'Color', 'blue', 'LineWidth', 2);
            grid on
            
%             title('Annual maxima','FontWeight','bold')
%             xlabel('Year [-]')
%             ylabel('SWE [mm]')
%             legend('observation', 'linear fit')
            title('Annual maxima','FontWeight','bold')
            xlabel('Year [-]')
            %ylabel('SWE [mm]')
            ylabel('s [kN/m^2]')
            legend('observation', 'linear fit')
            
            xlim([min(years)-1, max(years)+1])
            
        otherwise
            f_scatter = 0;
    end
    
end