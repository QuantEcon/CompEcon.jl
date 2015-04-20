function demapp10

%% DEMAPP10 Monopolist's Effective Supply Function

% Preliminary tasks
demosetup(mfilename)

%% Residual Function
resid = @(c,basis,p) ...
   p + funeval(c,basis,p)./(-3.5*p.^(-4.5)) ...
     - sqrt(funeval(c,basis,p)) ...
     - funeval(c,basis,p).^2;

% Approximation structure
n = 21; a = 0.5; b = 2.5;
basis = fundefn('cheb',n,a,b);
pnode = funnode(basis);

% Solve for effective supply function
c = [2;zeros(n-1,1)];
c = broyden(resid,c,basis,pnode);

% Setup plot
nplot = 1000;
pplot = nodeunif(nplot,a,b);
rplot = resid(c,basis,pplot);

% Plot effective supply
figure
plot(funeval(c,basis,pplot),pplot)
title('Monopolist''s Effective Supply Curve')
xlabel('Quantity')
ylabel('Price')

% Plot residual
figure
plot(pplot,rplot,pplot,0*pplot,'k--')  
title('Functional Equation Residual')
xlabel('Price')
ylabel('Residual')

% Save Plots as EPS Files
printfigures(mfilename,2)


% %% Residual Function
% function r = residual(c,basis,p)
% s = funeval(c,basis,p);
% k = sqrt(s)+s.^2;
% dprime = -3.5*p.^(-4.5);
% r = p + s./dprime - k;