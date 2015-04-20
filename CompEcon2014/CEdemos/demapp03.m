function demapp03

%% DEMAPP03 Basis functions and standard nodes for major approximation schemes

% Preliminary tasks
demosetup(mfilename)

% Set degree of approximation and endpoints of approximation interval
n = 9;
a = 0;
b = 1;

% Construct plotting grid
m = 1001;
x = nodeunif(m,a,b);

% Plot monomial basis functions
phi = zeros(m,n);
for j=1:n
    phi(:,j) = x.^(j-1);
end
basisplot(n,x,phi,'Monomial Basis Functions on [0,1].')

% Plot Chebychev basis functions and nodes
basis = fundefn('cheb',n,a,b);
phi   = funbase(basis,x);
xnode = funnode(basis);
basisplot(n,x,phi,'Chebychev Polynomial Basis Functions on [0,1].')
nodeplot(xnode,'Chebychev Nodes on [0,1]');

% Plot linear spline basis functions and nodes
basis = fundefn('spli',n,a,b,1);
phi   = full(funbase(basis,x));
xnode = funnode(basis);
basisplot(n,x,phi,'Linear Spline Basis Functions on [0,1].')
nodeplot(xnode,'Linear Spline Standard Nodes on [0,1].');

% Plot cubic spline basis functions and nodes
basis = fundefn('spli',n,a,b);
phi   = full(funbase(basis,x));
xnode = funnode(basis);
basisplot(n,x,phi,'Cubic Spline Basis Functions on [0,1].')
nodeplot(xnode,'Cubic Spline Standard Nodes on [0,1].');

% Save Plots as EPS Files
printfigures(mfilename,7)


%% Routine for plotting basis functions
function basisplot(n,x,phi,figtitle)
figure
fs = 18;
ymin = round(min(min(phi)));
ymax = round(max(max(phi)));
for j=1:n
    subplot(3,3,j); plot(x,phi(:,j),'LineWidth',5);
    axis([0 1 ymin-0.02 ymax+0.02]); 
    set(gca,'FontSize',fs,'box','off')
    set(gca,'Xtick',[0 1]);
    set(gca,'Ytick',[ymin ymax]);
    if j<7
        set(gca,'XtickLabel',[]);
    end
    if j~=1&j~=4&j~=7
        set(gca,'YtickLabel',[]);
    end
end
subplot(3,3,2); 
title(figtitle,'FontSize',14);


%% Routine for plotting approximation nodes
function nodeplot(xnode,figtitle)
figure
plot(xnode,zeros(size(xnode)),'*','LineWidth',5)
axis([-0.0001 1.0001 -0.05 0.05]); 
set(gca,'Xtick',[0 1]);
set(gca,'YTick',[],'DataAspectRatio',[4 2 2])
set(gca,'FontSize',18)
title(figtitle,'FontSize',14);