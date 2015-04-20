function demqua10

%% DEMQUA10 Monte Carlo Simulation of Time Series
%
% Simulate time series using Monte Carlo Method

% Preliminary tasks
demosetup(mfilename)


%% Perform Simulation
mu    = 0.005; 
sigma = 0.02;
n = 41;
m =  3;
z = randn(n,m);
logp = zeros(n,m);
logp(1,:) = log(2)*ones(1,m);
for i=1:n
  logp(i+1,:) = logp(i,:) + mu + sigma*z(i,:);
end


%% Plot Simulated Time Series
plot(exp(logp))
xlabel('Week')
ylabel('Price')
xlim([0 n])


%% Save Plots as EPS Files
printfigures(mfilename,1)