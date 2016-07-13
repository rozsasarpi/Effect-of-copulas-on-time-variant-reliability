% function mp_integral_test_02
clearvars
close all
clc

F = @(x,y) exp(-(x + y).^2);

theta = 2*pi;%pi/4;
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

Fuv = @(u,v) exp(-((1/sqrt(2)*u + -1/sqrt(2)*v) + (1/sqrt(2)*u + 1/sqrt(2)*v)).^2);

% function f = Fuv(u,v)
%     xy = [u,v]*R;
%     f = F(xy(1), xy(2));
% end

% masking
disp('-------------------------------------')
disp('zero masking dblquad')
tic
Q = dblquad(@(x,y)F(x,y).*(y <= x),-100,100,-100,100)
toc

disp('-------------------------------------')
disp('zero masking quad2d')
tic
Q = quad2d(@(x,y)F(x,y).*(y <= x),-100,100,-100,100,'MaxFunEvals',1000,'FailurePlot',true)
toc

% iterated
disp('-------------------------------------')
disp('home-made iterated')
tic
nn = 1000;
Iy = nan(nn,1);
xx = linspace(0,1,nn);
for ii = 1:nn
    Iy(ii) = integral(@(y) F(xx(ii), y) , -100, xx(ii));
end
plot(xx,Iy)
trapz(xx,Iy)
toc

disp('-------------------------------------')
disp('built-in iterated')
tic
Q = integral2(F,-100,100,-100,@(x)x, 'method', 'iterated')
toc

% rotated - rectangular
disp('-------------------------------------')
disp('rotated - rectangular')
tic
Q = dblquad(Fuv, -100,100, 0,-100)
toc

% end