%% DEMSCONT01 Ito Processes
%
% Simulate geometric Brownian motion, drift mu and volatility sigma

% Preliminary tasks
demosetup(mfilename)

% Model Parameters
T = 1;
n = 365;
h = T/n;
t = (0:h:T)';
mu = 0.1;
sigma = 0.05;

% Simulate
m = 3;
z = randn(n,m);
s = zeros(n+1,m);
s(1,:) = 1;
for i=1:n
  s(i+1,:) = s(i,:) + mu*s(i,:)*h + sigma*s(i,:)*sqrt(h).*z(i,:);
end

% Plot
figure
plot(t,s,'LineWidth',3)
xlabel('t')
ylabel('s(t)')
title('Simulated Ito Process, mu=0.1, sigma=0.05')

% Save Plots as EPS Files
printfigures(mfilename,1)