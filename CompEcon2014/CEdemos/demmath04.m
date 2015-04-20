function demmath04

%% DEMMATH04 Standard Copulas

% Preliminary tasks
demosetup(mfilename)

set(0,'defaulttextinterpreter','latex')

tau = 0.7;
nc = 5;
zlim = 2.5;
n = 800;
u = linspace(0,1,n);
[U1,U2] = meshgrid(u,u);
U = [U1(:) U2(:)];

figure
subplot(1,2,1)
hold on
c = copulapdf('Clayton',U,2*tau/(1-tau));
c = reshape(c,n,n);
z = icdf('Normal',u,0,1);
f = pdf('Normal',z,0,1);
c = c.*kron(f',f);
contour(z,z,c,nc,'Linewidth',2)
xlim([-zlim zlim]);
ylim([-zlim zlim]);
xlabel('$\tilde z_1$','FontSize',16)
ylabel('$\tilde z_2$','FontSize',16)
title('Clayton','FontSize',16)
xlabel('$z_1$')
ylabel('$z_2$')
axis square
subplot(1,2,2)
hold on
c = copulapdf('Gumbel',U,1/(1-tau));
c = reshape(c,n,n);
z = icdf('Normal',u,0,1);
f = pdf('Normal',z,0,1);
c = c.*kron(f',f);
contour(z,z,c,nc,'Linewidth',2)
xlim([-zlim zlim]);
ylim([-zlim zlim]);
xlabel('$\tilde z_1$','FontSize',16)
ylabel('$\tilde z_2$','FontSize',16)
title('Gumbel','FontSize',16)
xlabel('$z_1$')
ylabel('$z_2$')
axis square

figure
n = 100;
subplot(1,2,1)
hold on
u = copularnd('Clayton',2*tau/(1-tau),n);
z = u;
z = icdf('Normal',u,0,1);
scatter(z(:,1),z(:,2),'*')
xlim([-zlim zlim]);
ylim([-zlim zlim]);
xlabel('$\tilde z_1$','FontSize',16)
ylabel('$\tilde z_2$','FontSize',16)
title('Clayton','FontSize',16)
axis square
subplot(1,2,2)
hold on
u = copularnd('Gumbel',1/(1-tau),n);
z = u;
z = icdf('Normal',u,0,1);
scatter(z(:,1),z(:,2),'*')
xlim([-zlim zlim]);
ylim([-zlim zlim]);
xlabel('$\tilde z_1$','FontSize',16)
ylabel('$\tilde z_2$','FontSize',16)
title('Gumbel','FontSize',16)
axis square

% Save Plots as EPS Files
printfigures(mfilename,2)