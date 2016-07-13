clear all
close all
clc

%=================================================
%% OPTIONS

mu              = 1;                    % mean value
sigma           = 0.2*mu;               % standard deviation
distribution    = 'lnorm';              % distribution type, norm and lnorm are available

corr_length     = 6*1e5;                 % correlation length
pow             = 2;                    % exponent to calculate cov_mx

N               = 1000;                   % simulation number

A_mesh_x        = 0:20:800;             % node positions in x direction
A_mesh_y        = 0:200:6000;             % node positions in y direction
%=================================================

%=================================================
%% CALCULATION

%for lognormal distribution
mu_lognorm = log(mu.^2./sqrt((sigma).^2 + mu.^2));
sigma_lognorm = sqrt(log((sigma).^2./mu.^2 + 1));

[X, Y] = meshgrid(A_mesh_x, A_mesh_y);
Coord = [reshape(X,[],1), reshape(Y,[],1)];

n_x = length(A_mesh_x);
n_y = length(A_mesh_y);

cov_matrix = element_cov_matrix_2D(Coord, corr_length, pow);

% %Cholesky decomposition of cov_matrix
% L = chol(cov_matrix,'lower');
% % 
% ir1 = normrnd(mu, sigma, size(Coord,1), N);
% % creates dependent variables according to cov_matrix
% dr1 = L*ir1;
% % normalize the generated variables to get the desired mean of the marginal
% % pdfs
% dr1 = bsxfun(@rdivide, dr1, sum(L,2));
% 
% r = dr1;

% generates random numbers
switch lower(distribution)
    case {'normal', 'norm'}
        r = mvnrnd(ones(size(cov_matrix,1),1)*mu, cov_matrix*sigma^2, N)';
        
    case {'lognormal', 'lnorm'}
        r = MvLogNRand(ones(size(cov_matrix,1),1)*mu_lognorm, ones(size(cov_matrix,1),1)*sigma_lognorm, N, cov_matrix);
        r = r';
end

% aa = (diag(cov(r')))
% max(aa)
% min(aa)
% cov(r')/sigma
% mean(r,2)


%=================================================
%% VISUALIZATION 
%3D plot
figure('Units','Normalized','Position',[.15 .15 .7 .7])
subplot(2,3,1)

ZM = reshape(mean(r,2), size(X,1), size(X,2));
ZU = reshape(mean(r,2)+std(r,[],2), size(X,1), size(X,2));
ZB = reshape(mean(r,2)-std(r,[],2), size(X,1), size(X,2));

%............................................................
%one particular simulated surface - 3D
ZI = reshape(r(:,1), size(X,1), size(X,2));

subplot(2,3,1)
surf(X,Y,ZI)
%alpha 0.5
set(gca, 'DataAspectRatio', [repmat(max(diff(get(gca, 'XLim')), diff(get(gca, 'YLim'))), [1 2]) diff(get(gca, 'ZLim'))])
axis([0,A_mesh_x(end), 0,A_mesh_y(end), min(min(ZI)), max(max(ZI))])
% axis([0,A_mesh_x(end), 0,A_mesh_y(end), mu-2*sigma, mu+2*sigma])
title('Egy kiválasztott szimulált felület ábrája.')
xlabel('dimension x [mm]')
ylabel('dimension y [mm]')
zlabel('gap [mm]')

%............................................................
%mean +- sigma of simulated sufaces - 3D
subplot(2,3,2)
surf(X,Y,ZU)
hold on
surf(X,Y,ZB)
surf(X,Y,ZM)
alpha 0.5
set(gca, 'DataAspectRatio', [repmat(max(diff(get(gca, 'XLim')), diff(get(gca, 'YLim'))), [1 2]) diff(get(gca, 'ZLim'))])
axis([0,A_mesh_x(end), 0,A_mesh_y(end), mu-2*sigma, mu+2*sigma])
title('Szimulált felületek várható értéke és +- sigma konf. intervallumok.')
xlabel('dimension x [mm]')
ylabel('dimension y [mm]')
zlabel('gap [mm]')

%............................................................
%Contour plot
subplot(2,3,4)
%one particular simulated surface - contour
contourf(X,Y,ZI,10)
colorbar
axis equal
axis([0,A_mesh_x(end), 0,A_mesh_y(end)])
title('Egy kiválasztott szimulált felület ábrája (felsõvel azonos).')
xlabel('dimension x [mm]')
ylabel('dimension y [mm]')

%............................................................
% select the the first slice in y direction for 2D plot
r1 = r(1:size(X,1),:);
xx = Y(1:size(X,1));


for i=1:size(r1,2)
    rr(:,i) = interp1(xx,r1(:,i),linspace(xx(1),xx(end),size(r1,1)*3),'spline');
end

subplot(2,3,5)
%plot(rr','black','LineWidth',2)   %black
plot(rr(:, 1:min(size(r1,2),30)),'--','LineWidth',2)       %colored, maximum number of plots is limited!
grid on
hold on
plot(mean(rr,2),'blue','LineWidth',3)
plot(mean(rr,2)+std(rr,[],2),'red','LineWidth',3)
plot(mean(rr,2)-std(rr,[],2),'red','LineWidth',3)
axis([0 size(rr,1) 0.9*min(min(rr)) 1.1*max(max(rr))])
title('y irányban felvett metszetek (x=0 mellett) a szimulált felületekrõl.')
xlabel('dimension y [mm]')
ylabel('gap [mm]')

%............................................................
subplot(2,3,6)
%rrr = r(1,:);
rrr = reshape(r,[],1);
[counts,bins] = hist(rrr, min(size(r1,2),20)); % get counts and bin locations
barh(bins,counts)
axis([0 1.1*max(counts) 0.9*min(min(rr)) 1.1*max(max(rr))])
grid on
title('Gap hisztogramja a teljes felületre vonatkozóan.')
xlabel('gyakoriság [-]')
ylabel('gap [mm]')
%=================================================