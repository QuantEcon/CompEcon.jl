function demapp05

%% DEMAPP05 Chebychev polynomial and spline approximantion of various functions

% Preliminary tasks
demosetup(mfilename)

% Demonstrates Chebychev polynomial, cubic spline, and linear spline
% approximation for the following functions
%   1: y = sqrt(x+1);
%   2: y = exp(-x);
%   3: y = 1./(1+25*x.^2).
%   4: y = sqrt(abs(x));

% Set degree of approximation and endpoints of approximation interval
n =  7;                       % degree of approximation
a = -1;                       % left endpoint
b =  1;                       % right endpoint

% Construct Chebychev, cubic spline, and linear spline bases
chebbasis = fundefn('cheb',n,a,b);
csplbasis = fundefn('spli',n,a,b);
lsplbasis = fundefn('spli',n,a,b,1);

% Construct uniform grid for error ploting
x = nodeunif(2001,a,b);

for ifunc=1:4
    
   % Construct interpolants
   ccheb = funfitf(chebbasis,@f,ifunc);
   ccspl = funfitf(csplbasis,@f,ifunc);
   clspl = funfitf(lsplbasis,@f,ifunc);
   
   % Compute actual and fitted values on grid
   y = f(x,ifunc);                          % actual
   ycheb = funeval(ccheb,chebbasis,x);      % Chebychev approximant
   ycspl = funeval(ccspl,csplbasis,x);      % cubic spline approximant
   ylspl = funeval(clspl,lsplbasis,x);      % linear spline approximant

   % Plot function approximations
   figure
   ymin = floor(min(y)); ymax=ceil(max(y)); 
   if ifunc==3
     axislim=[a b -0.2 1.2];
   else
     axislim=[a b ymin ymax];
   end
   subplot(2,2,1); plot(x,y,'k');     axis(axislim); title('Function');
   subplot(2,2,2); plot(x,ycheb,'k'); axis(axislim); title('Chebychev');
   subplot(2,2,3); plot(x,ycspl,'k'); axis(axislim); title('Cubic Spline');
   subplot(2,2,4); plot(x,ylspl,'k'); axis(axislim); title('Linear Spline'); 
   
   % Plot function approximation error
   figure
   ymin=floor(min(y)); ymax=ceil(max(y)); 
   if ifunc==3
     axislim=[a b -0.2 1.2];
   else
     axislim=[a b ymin ymax];
   end
   plot(x,[ycheb-y ycspl-y],'LineWidth',4)
   xlabel('x')
   ylabel('Approximation Error')
   legend('Chebychev Polynomial','Cubic Spline')  
   legend('boxoff')
   
end

% Save Plots as EPS Files
printfigures(mfilename,8)


%% Functions to be approximated
function y = f(x,i)
switch i
   case 1, y = 1 + x + 2*x.^2 - 3*x.^3;
   case 2, y = exp(-x);
   case 3, y = 1./(1+25*x.^2);
   case 4, y = sqrt(abs(x));
end