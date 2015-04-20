function demslv09

%% DEMSLV09 Hard nonlinear complementarity problem with Billup's function
%
% Solve hard nonlinear complementarity problem on R using semismooth and 
% minmax methods.  Problem involves Billup's function.  Minmax formulation 
% fails; semismooth formulation suceeds.

% Preliminary tasks
demosetup(mfilename)
warning off

% Figure parameters
fs = 12;  % FontSize
bs = 18;  % BulletSize

% Generate problem test data
xinit = 0;
a = 0;
b = inf;

% Print table header
fprintf('Hundreds of seconds required to solve hard nonlinear complementarity\n')
fprintf('problem using Newton semismooth and minmax formulations\n')
fprintf('Algorithm           Time      Norm         x\n');

% Solve by applying Newton method to minmax formulation
tic
optset('ncpsolve','type','minmax');
[x,z] = ncpsolve(@billups,a,b,xinit);
fprintf('Newton minmax     %6.2f  %8.0e     %5.2f\n',100*toc,norm(minmax(x,a,b,z)),x)

% Solve by applying Newton method to semismooth formulation
tic
optset('ncpsolve','type','ssmooth');
[x,z] = ncpsolve(@billups,a,b,xinit);
fprintf('Newton semismooth %6.2f  %8.0e     %5.2f\n',100*toc,norm(minmax(x,a,b,z)),x)

figure

subplot(1,2,1)
hold on
x = nodeunif(500,-0.5,2.5);
plot(x,billupss(x),x,billupsm(x))
legend('Semismooth','Minmax','Location','S')
legend('boxoff')
plot(x,zeros(size(x)),'k:');
set(gca,'FontSize',fs)
xlabel('x')   
axis([-0.5 2.5 -1.0 1.5],'square')
title('Difficult NCP')

subplot(1,2,2)
hold on
x = nodeunif(500,-0.03,0.03);
plot(x,billupss(x),x,billupsm(x))
legend('Semismooth','Minmax','Location','S')
legend('boxoff')
plot(x,zeros(size(x)),'k:');
set(gca,'FontSize',fs)
xlabel('x') 
axis([-.03 .03 -.04 .06],'square')
title('Difficult NCP Magnified')

% Save Plots as EPS Files
printfigures(mfilename,1,1,0)


%% Billups' function
function [fval,fjac] = billups(x)
fval = 1.01-(1-x).^2;
fjac = 2*(1-x);


%% Minmax reformulation of Billups' function
function [fval,fjac] = billupsm(x)
[fval,fjac] = billups(x);
if nargout<2
   fval = minmax(x,0,inf,fval,fjac);
else
   [fval,fjac] = minmax(x,0,inf,fval,fjac);
end


%% Semismooth reformulation of Billups' function
function [fval,fjac] = billupss(x)
[fval,fjac] = billups(x);
if nargout<2
   fval = ssmooth(x,0,inf,fval,fjac);
else
   [fval,fjac] = ssmooth(x,0,inf,fval,fjac);
end