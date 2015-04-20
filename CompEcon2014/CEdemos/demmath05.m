function demmath05

%% DEMMATH05 Illustrates Implicit Function Theorem

% Preliminary tasks
demosetup(mfilename)


%% Implicit Function Theorem

n = 100;
x = nodeunif(n,0,2);
F = @(y) y.^5+y.^3-x.^2-1;
y = ones(n,1);
y = broyden(F,y);
f1 = 1/4;
f2 = 3/64;
y1 = 1+f1*(x-1);
y2 = y1 + 0.5*f2*(x-1).^2;

figure
plot(x,y,x,y1,x,y2)
title('Illustration of Implicit Function Theorem')
legend('Implicit Function','1st Order Approximation','2nd Order Approximation','Location','Best')
xlabel('x')
ylabel('y')
xlim([0 2])
set(gca,'xTick',[0 1 2])
ylim([0.7 1.3])
set(gca,'yTick',[0.7 1.0 1.3])


%% Inverse Function Theorem

n = 200;
p = nodeunif(n,0.4,3.0);
q = 0.5*p.^-0.2 + 0.5*p.^-0.5;
g1 = -1/0.35;
g2 = 0.495/0.35^3;
p1 = 1+g1*(q-1);
p2 = p1 + 0.5*g2*(q-1).^2;
figure
plot(q,p,q,p1,q,p2)
title('Illustration of Inverse Function Theorem')
legend('Invese Function','1st Order Approximation','2nd Order Approximation','Location','Best')
xlabel('q')
ylabel('p')
xlim([0.7 1.3])
set(gca,'xTick',[0.7 1.0 1.3])
set(gca,'yTick',[0 1 2 3])

%% Save Plots as EPS Files
printfigures(mfilename,2)