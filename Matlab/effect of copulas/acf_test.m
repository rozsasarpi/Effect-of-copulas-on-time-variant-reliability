clear all
close all
clc

rng(1); % For reproducibility
Mdl = arima('MA',{-0.5 0.4},'Constant',0,'Variance',1);

y = simulate(Mdl,1000);

[ACF,lags,bounds] = autocorr(y,[],2);
bounds

autocorr(y)