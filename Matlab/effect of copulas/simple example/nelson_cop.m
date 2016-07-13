function nelson_cop

theta = 10;

%Nelson #13
I = @(t) inelson(t);

x = 0:0.01:1;
y = I(x);
plot(x,y)

%Clayton
% I = @(t) 1/theta*(t.^-theta-1)./(-t.^(-theta-1));

IN = integral(I, 0,1, 'AbsTol', 1e-12, 'RelTol', 1e-6);

tau = 1 + 4*IN

%Clayton - closed form solutioon
theta/(2+theta)


function I = inelson(t) 

I = nan(size(t));
idx = t == 0;
I(~idx) = ((1-log(t(~idx))).^theta-1)./(-(theta.*(1-log(t(~idx))).^(theta-1)))./t(~idx);
I(idx) = -Inf;

end

end