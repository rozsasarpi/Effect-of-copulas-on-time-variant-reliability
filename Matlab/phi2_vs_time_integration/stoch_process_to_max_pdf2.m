clear all
close all
clc

% OPTIONS

% number of realizations per block (snow events per year, independent)
nr = 25;

%number of samples, blocks (years with observation)
nb = 1e6;

mu = 3500;
sigma = 700;

% SIMULATION
R = normrnd(mu, sigma, nr, nb);

maxR = max(R,[],1);

% PLOT
[f,x] = hist(maxR, 50);
bar(x,f/trapz(x,f))
hold on
%xlim([min(x)*0.95, max(x)*1.05])

xx = min(x)*0.95:0.1:max(x)*1.05;
param= evfit(-maxR);
[M, V] = evstat(-param(1), param(2));
M
sqrt(V)
p = evpdf(-xx,param(1),param(2));
plot(xx,p,'color','r')

meanX = mean(maxR);
stdX = std(maxR);
p2 = gumbelpdf(xx, meanX, stdX);
plot(xx,p2,'color','g')

legend(['max(X_1,X_2,..X_{',num2str(nb),'})'], 'MLE-Gumbel', 'MOM-Gumbel')