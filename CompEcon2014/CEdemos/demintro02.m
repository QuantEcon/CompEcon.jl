function demintro02

%% DEMINTRO02 Rational Expectations Agricultural Market Model

% Preliminary tasks
demosetup(mfilename)

% Generate yield distribution
sigma = 0.2;
[y,w] = qnwlogn(25,-0.5*sigma^2,sigma^2);

% Compute rational expectations equilibrium using function iteration,
% iterating on acreage planted
ptarg = 1;
a = 1;
for it=1:50
  aold = a;
  a = 0.5 + 0.5*w'*max(1.5-0.5*a*y,ptarg);
  fprintf('%3i %8.4f %8.1e\n',[it a norm(a-aold)])
  if norm(a-aold)<1.e-8, break, end
end

% Intermediate outputs
q = a*y;            % quantity produced in eaxxch state
p = 1.5-0.5*a*y;    % market price in each state
f = max(p,ptarg);   % farm price in each state
r = f.*q;           % farm revenue in each state
g = (f-p).*q;       % government expenditures
[xavg,xstd] = discmoments(w,[p f r g])

% Print results 
fprintf('Variable                   Expect     Std Dev\n')
fprintf('Market Price             %8.4f   %8.4f\n',xavg(1),xstd(1))
fprintf('Farm Price               %8.4f   %8.4f\n',xavg(2),xstd(2))
fprintf('Farm Revenue             %8.4f   %8.4f\n',xavg(3),xstd(3))
fprintf('Government Expenditures  %8.4f   %8.4f\n',xavg(4),xstd(4))

% Generate fixed-point mapping
aeq = a;
a = nodeunif(100,0,2);
g = zeros(100,1);
for i=1:100
  g(i) = 0.5 + 0.5*w'*max(1.5-0.5*a(i)*y,1);
end

% Graph rational expectations equilibrium
figure
hold on
plot(a,g,'LineWidth',5)
plot(a,a,'k--','LineWidth',3)
plot([0   aeq],[aeq aeq],'r--','LineWidth',3)
plot([aeq aeq],[0   aeq],'r--','LineWidth',3)
set(gca,'xtick',[0 1 2])
set(gca,'ytick',[0 1 2])
axis([0 2 0 2],'square')
xlabel('Acreage Planted')
ylabel('Rational Acreage Planted')
ftext(0.05,1.35,'g(a)','left','middle')
ftext(0.25,0.00,'45^o','center','bottom')
set(gca,'XTick',[0 aeq 2])
set(gca,'YTick',[0 aeq 2])
set(gca,'XTickLabel',{'0' 'a*' '2'})
set(gca,'YTickLabel',{'0' 'a*' '2'})
bullet(aeq,aeq,24)

% % Compute rational expectations equilibrium using function iteration,
% % iterating on expected farm price
% Ef = 1;
% for it=1:50
%   Efold = Ef;
%   Ef = w'*(max(1.5-0.5*(0.5 + 0.5*Ef)*y,ptarg));
%   if abs(Ef-Efold)<1.e-10, break, end
% end

% Compute rational expectations equilibrium as a function of the target
% price
nplot = 50;
ptarg = nodeunif(nplot,0,2);
a = 1;
for ip=1:nplot
  for it=1:50
    aold = a;
    a = 0.5 + 0.5*w'*max(1.5-0.5*a*y,ptarg(ip));
    if norm(a-aold)<1.e-10, break, end
  end
  q = a*y;              % quantity produced
  p = 1.5-0.5*a*y;      % market price
  f = max(p,ptarg(ip)); % farm price
  r = f.*q;             % farm revenue
  g = (f-p).*q;         % government expenditures
  [xavg,xstd] = discmoments(w,[p f r g]);
  Ep(ip) = xavg(1);
  Ef(ip) = xavg(2);  
  Er(ip) = xavg(3);
  Eg(ip) = xavg(4);
  Sp(ip) = xstd(1);
  Sf(ip) = xstd(2);  
  Sr(ip) = xstd(3);
  Sg(ip) = xstd(4);
end

% Graph expected prices vs target price
figure
hold on
plot(ptarg,Ep,ptarg,Ef,'LineWidth',5)
plot(ptarg([1 nplot]),Ep([1 1]),'k:','LineWidth',3)
title('Expected Prices vs Target Price')
xlabel('Target Price')
ylabel('Expectation')
legend('Market Price','Farm Price','Location','NW')
legend boxoff
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0.5 1.0 1.5 2.0])
set(gca,'YTickLabel',{'0.5' '1.0' '1.5' '2.0'})
ylim([0.5 2])

% Graph price standard deviations vs target price
figure
hold on
plot(ptarg,Sp,ptarg,Sf,'LineWidth',5)
plot(ptarg([1 nplot]),Sf([1 1]),'k:','LineWidth',3)
title('Price Variabilities vs Target Price')
xlabel('Target Price')
ylabel('Standard Deviation')
legend('Market Price','Farm Price','Location','NW')
legend boxoff
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 0.01 0.02])

% Graph expected farm revenue vs target price
figure
hold on
plot(ptarg,Er,'LineWidth',5)
plot(ptarg([1 nplot]),Er([1 1]),'k:','LineWidth',3)
title('Expected Farm Revenue vs Target Price')
xlabel('Target Price')
ylabel('Expectation')
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 1 2 3])
set(gca,'YTickLabel',{'0' '1' '2' '3'})
ylim([0.8 3])

% Graph standard deviation of farm revenue vs target price
figure
hold on
plot(ptarg,Sr,'LineWidth',5)
plot(ptarg([1 nplot]),Sr([1 1]),'k:','LineWidth',3)
title('Farm Revenue Variability vs Target Price')
xlabel('Target Price')
ylabel('Standard Deviation')
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 0.2 0.4])
set(gca,'YTickLabel',{'0' '0.2' '0.4'})

% Graph expected government expenditures vs target price
figure
hold on
plot(ptarg,Eg,'LineWidth',5)
plot(ptarg([1 nplot]),Eg([1 1]),'k:','LineWidth',3)
title('Expected Government Expenditures vs Target Price')
xlabel('Target Price')
ylabel('Expectation')
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 1 2])
ylim([-0.05 2])


%% Save Plots as EPS Files
printfigures(mfilename,6)