clear all
close all
clc

N = 1e4;
pf = 0.01;
I = rand(N,1)<pf;

nn = 1:N;
sd = zeros(1,N);

for i = 1:N
    sd(i) = std(I(1:i));
end

disp('probability of failure (final):')
mean(I)

se = sd.^2./(nn).^(1/2);
plot(nn,sd)
xlabel('simulation number')
ylabel('std of indicator rv')
disp('standard deviation of indicator rv (final):')
sd(end)
disp('coeff. of var. of indicator rv (final):')
sd(end)/mean(I)
%aa = (pf^2*(1-pf) + pf*(1-pf)^2)/pf %1-pf

%% SE
figure
plot(nn,se)
xlabel('simulation number')
ylabel('standard error of failure estimate')
disp('standard error of failure estimate (final):')
se(end)

%% COV
figure
plot(nn,se)
xlabel('simulation number')
ylabel('coeff. of var. of failure estimate')
disp('coeff. of var. of failure estimate (final):')
se(end)/mean(I)
(1-pf)/sqrt(N)


