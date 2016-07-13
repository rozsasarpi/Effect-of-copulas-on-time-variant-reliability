% simple random walk test

% dimensions and discretization of state space
x = 1:10;
y = 1:10;

% construction of transition matrix

%% most simple random walk without transition matrix
x0      = [0, 0];

for i = 1:N
% random direction
alpha   = 2*pi*rand(numel(x0)-1, 1);

% random or discrete step size
step    = 1;

% new position
x(i,1)  = x(i-1,1) + cos(alpha)*step;
x(i,2)  = x(i-1,2) + cos(alpha)*step;

