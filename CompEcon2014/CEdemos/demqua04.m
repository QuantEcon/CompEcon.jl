function demqua04

%% DEMQUA04 Area under normal pdf using Simpson's rule

% Preliminary tasks
demosetup(mfilename)

n = 11; a = 0; z = 1;
[x,w] = qnwsimp(n,a,z);
prob = 0.5 + w'*f(x)

figure
hold on
b = 4;
a = -b;
n = 500;
x = nodeunif(n,a,b);
y = f(x);
for i=1:n
    if x(i)<z, plot ([x(i) x(i)],[0 y(i)],'y','LineWidth',2), end
end
plot([a b],[0 0],'k-')
plot([z z],[0 f(z)],'k-','LineWidth',2)
plot(x,y)
set(0,'defaulttextinterpreter','latex')
ftext( 1,-0.01,'$z$','center','top',28);
ftext(-2, 0.20,'$\Pr\left(\tilde Z\leq z\right)$','right','middle',28);
annotation('textarrow',[0.33 0.45],[0.52 0.36]);
set(gca,'ytick',[0.0 0.2 0.4])
axis off

% Save Plots as EPS Files
printfigures(mfilename,1)


function y=f(x)
y = sqrt(1/(2*pi))*exp(-0.5*x.^2);