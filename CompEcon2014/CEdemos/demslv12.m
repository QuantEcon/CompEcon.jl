function demapp12

%% DEMAPP12  Convergence rates for different NLP methods

% Preliminary tasks
demosetup(mfilename)

n = 30;
xinit = 2;

x = xinit;
err1 = [];
x1 = [];
for it=1:n
    x1 = [x1;x];
    err1 = [err1; log10(x)];
    f = exp(x)-1;
    J = exp(x);
    s = -f/J;
    x = x + s;
end

x    = xinit;
xold = xinit+sqrt(eps);
fold = exp(xold)-1;
x2 = [];
err2 = [];
for it=1:n
    x2 = [x2;x];
    err2 = [err2; log10(x)];
    f = exp(x)-1;
    s = -f/((f-fold)/(x-xold));
    xold = x;
    fold = f;
    x = x + s;
end

x = xinit;
x3 = [];
err3 = [];
for it=1:n
    x3 = [x3;x];
    err3 = [err3; log10(abs(x))];
    x = x - exp(x) + 1;
end

figure
hold on
plot(1:length(err1),err1,1:length(err2),err2,1:length(err3),err3)
set(gca,'box','off')
xlabel('Iteration')
ylabel('Log10 Error')
legend('Newton''s Method','Broyden''s Method','Function Iteration','Location','SouthWest')
xlim([0 12])
ylim([-15 2])
legend('boxoff')

% Save Plots as EPS Files
printfigures(mfilename,1)