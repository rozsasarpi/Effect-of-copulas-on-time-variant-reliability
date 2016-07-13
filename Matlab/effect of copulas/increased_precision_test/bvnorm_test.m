clearvars
close all
clc

rho     = 1-1e-5;
idig    = 50;

x       = [10,10];
mpx     = mp(x);

P1      = 1-bivnor(-mpx(1),-mpx(2),rho,idig)
P2      = 1-bivnor(-mpx(1),-mpx(2),rho)
% P3      = 1-mvncdf([x(1),x(2)],[0,0],[1, rho; rho, 1])
P3      = 1-mvncdf([x(1),x(2)],[0,0],[1, rho; rho, 1])

(P1-P2)/P1

%%
u = mp('[0.5, 0.5]');
C1 = binorm_copulacdf_mp(u, rho)
C2 = copulacdf('Gaussian', double(u), rho)

(C1-C2)/C1

% TEST PLOT
% nn = 20; mm = nn;
% 
% xx = linspace(-3,3,nn);
% yy = linspace(-3,3,mm);
% P  = nan(nn,mm);
% P2 = P;
% for ii = 1:nn
%     for jj = 1:mm
%         P(ii,jj) = bivnor(-xx(ii),-yy(jj),rho);
%         P2(ii,jj)= mvncdf([xx(ii),yy(jj)],[0,0],[1, rho; rho, 1]);
%     end
% end
% 
% [X,Y] = meshgrid(xx,yy);
% surfc(X,Y,P)
% figure
% surfc(X,Y,P2)
