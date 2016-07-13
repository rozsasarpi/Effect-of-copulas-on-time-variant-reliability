clearvars
close all
clc

mp.Digits(50);

u = mp('[0.99, 0.99]');
theta = mp('100000');
copulacdf('Clayton', double(u), double(theta))

(u(1)^-theta + u(2)^-theta-1)^(-1/theta)

biclay_copulacdf_mp(u, theta)