function demqua05

%% DEMQUA05 Willingness to pay, expected utility model

% Preliminary tasks
demosetup(mfilename)

n = 100;
mu = 0;     
var = 0.1;
alpha = 2;  
ystar = 1;
[y,w] = qnwlogn(n,mu,var);
expectedutility = -w'*exp(-alpha*y)
certainutility  = -exp(-alpha*ystar)

ystar = -log(-expectedutility)/alpha;
wtp = w'*y-ystar