% Probability distribution function (pdf) of the product of two independent random variables
%
% Z = X1*X2
%
% pdz = PROD_2RV(pd1, pd2)

function pdz = prod_2rv(pd1, pd2)

fz = @(z) integral(@(x) pd1.fx_fun(x).*pd2.fx_fun(z./x).*1./abs(x), 0, inf);

pdz.fx_fun = fz;

end


% pd1.fx_fun = @(x) lognpdf(x,1,0.1);
% pd2.fx_fun = @(x) lognpdf(x,1,0.1);

% x = 0:0.1:20;
% 
% tic
% y = zeros(numel(x),1);
% for i = 1:numel(x)
%     y(i) = fz(x(i));
% end
% toc
% 
% plot(x,fx1_fun(x))
% hold on
% plot(x,y,'r')