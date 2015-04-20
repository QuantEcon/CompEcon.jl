function demintro01

%% DEMINTRO01 Inverse Demand Problem

% Preliminary tasks
demosetup(mfilename)

% Compute price that will clear market of q=2 using Newton's method
p = 0.25;
for it=1:100
    f = 0.5*p^-0.2 + 0.5*p^-0.5 - 2;
    d = -0.01*p^-1.2 - 0.25*p^-1.5;
    s = -f/d;
    p = p + s;    
    fprintf('iteration %3i  price %8.4f\n',[it p])
    if norm(s)<1.e-8, break, end
end

% Generate demand function
pstar = p;
qstar = 0.5*pstar^-0.2 + 0.5*pstar^-0.5;
n = 100;
a = 0.02;
b = 0.40;
p = nodeunif(n,a,b);
q = 0.5*p.^-0.2 + 0.5*p.^-0.5;

% Graph demand and inverse demand functions
figure
subplot(1,2,1)
plot(p,q)
set(gca,'xtick',[0.0 0.2 0.4])
set(gca,'ytick',[0 2 4])
set(gca,'box','off')
axis([0 b 0 4],'square')
title('Demand')
xlabel('p')
ylabel('q')
subplot(1,2,2)
hold on
plot(q,p)
plot([    0 qstar],[pstar pstar],'r--','LineWidth',3)
plot([qstar qstar],[0     pstar],'r--','LineWidth',3)
ftext(0,pstar,'p^*','right','middle');
bullet(2,pstar,24)
set(gca,'ytick',[0.0 0.2 0.4])
set(gca,'xtick',[0 2 4])
set(gca,'box','off')
axis([0 4 0 b],'square')
title('Inverse Demand')
xlabel('q')
ylabel('p')

% Save Plots as EPS Files
printfigures(mfilename,1,1,0)