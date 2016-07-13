clearvars
close all
clc
format long

k_tau       = 0.8;

theta       = 2*k_tau/(1-k_tau);
clay_pdf    = @(u1, u2) reshape(copulapdf('Clayton', [reshape(u1,[],1), reshape(u2,[],1)], theta), size(u1,1), []);
cclay_pdf   = @(u1, u2) reshape(copulapdf('Clayton', [reshape(1-u1,[],1), reshape(1-u2,[],1)], theta), size(u1,1), []);

U = [0.8, 0.8];

% Clayton
disp('========== Clayton ==========')
copulacdf('Clayton', U, theta)

integral2(clay_pdf, 0,U(1), 0,U(2))


% rotClayton
disp('========== rotClayton ==========')
U(1) + U(2) - 1 + copulacdf('Clayton', 1-U, theta)

integral2(cclay_pdf, 0,U(1), 0,U(2))


