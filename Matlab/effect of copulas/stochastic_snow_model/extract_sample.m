clearvars; close all; clc

load('D:\Working folder\Matlab working folder\snow load\obs_database\obs_database_swe.mat')
%database_code: lon(E), lat(N), country(1-hun), altitude
load('d:\Working folder\Matlab working folder\snow load\obs_database\database_code.mat','database_code')

plot_obs    = 1;
%-----------------
E = 19.1;
N = 47.5;
% E = 20.6;
% N = 47.8;
%-----------------

idx_E       = database_code(:,1) == E;
idx_N       = database_code(:,2) == N;
idx         = logical(idx_E.*idx_N);
point_ID    = find(idx == 1);
% point_ID    = 2547; % Budapest

x           = obs_database(2:end,point_ID);
i           = find(x==9999, 1, 'last' )+1;
% x           = x(i:end);
x(1:i)      = NaN;
x(x == -1)  = NaN;

cmp         = get(groot,'defaultAxesColorOrder');

%%
close all
% clear all
clc
% i           = 14;
% idx1        = (1+(i-1)*365);
% idx2        = (200+(i-1)*365);
% % plot(x(idx1:idx2))
% % plot(x(1:10*365))
%
% % E = 20.6; N = 47.8;
% nn = ceil(length(x)/7);
%
% t = 1:length(x);
% % t = idx1:idx2;
% s = x(t);

s = x;

[pks,locs] = findpeaks(s);

% remove the smaller peak if they are too close (th close)
locs_r  = locs;
pks_r   = pks;
flag    = 0;
th      = 4; % days

counter = 0;
while any(diff(locs_r) <= th)
    idx     = find(diff(locs_r) <= th);
    [~,I]   = min([pks_r(idx), pks_r(idx+1)],[],2);
    
    pks_r(idx+I-1)  = [];
    locs_r(idx+I-1) = [];
    counter = counter + 1;
end
disp(['counter: ', num2str(counter)])

% plot(t, s)
% plot(locs_r, pks_r)
% hold on
yl      = 0;
yu      = 80;
n_year  = 49; % number of complete winter seasons
np      = nan(n_year,1);
sample  = cell(n_year,1);
X       = cell(n_year,1);
yy      = 0:365.25:length(s);
figure('Position', [100, 600, 600, 200])
for ii = 1:n_year
    idx1        = floor(365.25/2+(ii-1)*365.25);
    idx2        = ceil(idx1 + 365.25);
    
    idx         = locs_r > idx1 & locs_r < idx2;
    if mod(ii,2) == 1
        area([idx1, idx1, idx2, idx2], [yl, yu, yu, yl],...
            'FaceColor', 0.9*[1,1,1],...
            'EdgeColor', 'none')
    end
    hold on
    plot(locs_r(idx), pks_r(idx),'o-', 'Color', cmp(2,:))
    h = plot(idx1:idx2, s(idx1:idx2), 'black');
    
    xl = max(idx1-4*365.25,0);
    xu = idx2;
    xlim([xl, xu])
    ylim([yl, yu])
    xt = yy(yy > xl & yy < xu);
    set(gca, 'XTick', xt)
    set(gca, 'XTickLabel', 1961 + xt/365.25)
    ylabel('SWE [mm]')
    xlabel('Year [-]')
    box off
    plot([xl, xl], ylim, 'w', 'Linewidth', 2)
    set(gca,'TickLabelInterpreter', 'LaTeX')

    np(ii)      = sum(idx);
    sample{ii}  = pks_r(idx);
    X{ii}       = locs_r(idx);
    
    if ii == 22
%         ax1 = gca;
%         set(ax1,'box','off')
%         ax2 = axes('Position', get(ax1, 'Position'),'Color','none');
%         set(ax2,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')
        
%         keyboard
    end
end

idx = np < 1;
np(idx) = [];
sample(idx) = [];
X(idx) = [];

% figure
% hist(np)

save('snow_ensembles.mat', 'sample', 'X')
