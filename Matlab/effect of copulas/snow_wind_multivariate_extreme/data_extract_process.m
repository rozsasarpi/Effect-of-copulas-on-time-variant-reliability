
% CONSIDER MORE THAN ONE DAY!! SNOW-WIND

clear all
close all
clc
%Budapest [18.9, 19.3], [47.6, 47.3] -> 12 grid points
%idx;   Elong;      Nlat
%2546   19          47.5

% for exploration
% Elong   = 18.9:0.1:19.3;
% Nlat    = 47.3:0.1:47.6;

% seleceted grid point
Elong = 18.9;
Nlat  = 47.5;

[E, N] = meshgrid(Elong, Nlat);

Folder = 'D:\Working Folder\Matlab working folder\snow load\obs_database';
load([Folder,'\database_code.mat'],'database_code')

for i = 1:numel(E)
    
iE = database_code(:,1) == E(i);
iN = database_code(:,2) == N(i);
station_idx = find(iE.*iN == 1);

fig = 'n';      % plot figures?

%====================================
%
%====================================


years = (1960:2010)';

% snow (-9999 replaced with 1e-3)
%....................... should be loaded only once
load([Folder,'\obs_database.mat'],'obs_database')
% load([Folder,'\database_code.mat'],'database_code')
obs_database_swe = obs_database;
%.......................

rowNumber = size(obs_database,1);
obs_data_swe  = obs_database_swe(:,station_idx);
[annual_max_swe, annual_max_swe_idx, ~] = annual_max(obs_data_swe, rowNumber, years, fig);


% wind (-9999 kept)
%....................... should be loaded only once
load([Folder,'\obs_database_ws10.mat'],'obs_database')
% load([Folder,'\database_code_ws10.mat'],'database_code')
obs_database_ws10 = obs_database;
%.......................

rowNumber       = size(obs_database,1);
obs_data_ws10   = obs_database_ws10(:,station_idx);
% [annual_max_ws10, ~] = annual_max(obs_data_ws10, rowNumber, years, fig);
accomp_ws10     = obs_data_ws10(annual_max_swe_idx);

% remove 1960 with limited observations
annual_max_swe  = annual_max_swe(2:end);
accomp_ws10     = accomp_ws10(2:end);
years           = years(2:end);

% dependence measures
ro      = corr(annual_max_swe, accomp_ws10);
k_tau   = 2*asin(ro)/pi;

disp(station_idx)
disp(['Pearson correlation coefficient (ro): ', num2str(ro)])
disp(['Kendall rank correlation (tau): ', num2str(k_tau)])
disp('=============')

plot(years, annual_max_swe/max(annual_max_swe), '-o')
hold on
plot(years, accomp_ws10/max(accomp_ws10), 'r-o')

% pseudo-observations
n = numel(annual_max_swe);
u_swe   = tiedrank(annual_max_swe)/(n+1);
u_ws10  = tiedrank(accomp_ws10)/(n+1);
figure
plot(u_swe, u_ws10,'o')
xlabel('u_{SWE}')
ylabel('u_{WS10}')


end