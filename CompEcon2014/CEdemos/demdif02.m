function demdif02

%% DEMDIF02
% Demonstrates accuracy of one- and two-sided finite-difference derivatives 
% of e(x) at x=1 as a function of step size h.

% Preliminary tasks
demosetup(mfilename)

n = 21;
c = -nodeunif(n,0,15);
c = sort(c)
h = 10.^c;
x = 1;

% One-sided finite difference derivative
u = x+h;
l = x;
d1 = (exp(u)-exp(l))./(u-l)-exp(x);
d1 = log10(abs(d1));
% d1 = smooth(c,d1,0.2,'rloess');
e1 = log10(eps^(1/2));

% Two-sided finite difference derivative
u = x+h;
l = x-h;
d2 = (exp(u)-exp(l))./(u-l)-exp(x);
d2 = log10(abs(d2));
% d2 = smooth(c,d2,0.2,'rloess');
e2 = log10(eps^(1/3));

% Plot finite difference derivatives
figure
hold on
plot(c,d1,c,d2)
plot([e1 e1],[-15 5],':','linewidth',2)
plot([e2 e2],[-15 5],':','linewidth',2)
xlim([-15 0])
xlabel('Log_{10}(h)')
ylabel('Log_{10} Approximation Error')
set(0,'defaulttextinterpreter','latex')
text(e1,2,'$\sqrt{\epsilon}$','fontsize',18)
text(e2,2,'$\sqrt[3]{\epsilon}$','fontsize',18)
legend('One-Sided','Two-Sided','Location','SW')
title('Error in Numerical Derivatives')

% Save Plots as EPS Files
printfigures(mfilename,1)