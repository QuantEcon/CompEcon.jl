function demapp04

%% DEMAPP04 
% Uniform-node and Chebychev-node polynomial approximation of Runge's function
% and compute condition numbers of associated interpolation matrices.

% Preliminary tasks
demosetup(mfilename)
warning off

% Set endpoints of approximation interval
a =  -1;                            % left endpoint
b =   1;                            % right endpoint

% Construct plotting grid
nplot = 1001;
x = nodeunif(nplot,a,b);
y = runge(x);

% Plot Runge's Function
figure
plot(x,y)
xlabel('x'); ylabel('y');
title('Runge''s Function');

% Initialize data matrices
n=3:2:40;
errunif = zeros(nplot,length(n));
errcheb = zeros(nplot,length(n));
nrmunif = zeros(length(n),1);
conunif = zeros(length(n),1);
nrmcheb = zeros(length(n),1);
concheb = zeros(length(n),1);

% Compute approximation errors on refined grid and 
% interpolation matrix condition numbers
for i=1:length(n)
   % Uniform-node monomial-basis approximant
   xnodes = nodeunif(n(i),a,b);  
   c = polyfit(xnodes,runge(xnodes),n(i)-1);
   yfit = polyval(c,x);
   phi = ones(n(i),1);      
   for j=1:n(i)-1
      phi = [phi xnodes.^j];
   end
   errunif(:,i) = yfit-y;
   nrmunif(i) = log10(norm(yfit-y,inf));
   conunif(i) = log10(cond(phi,2));
   % Chebychev-node Chebychev-basis approximant
   basis = fundefn('cheb',n(i),a,b);  
   c = funfitf(basis,@runge);          
   yfit = funeval(c,basis,x);   
   phi = funbase(basis);   
   errcheb(:,i) = yfit-y;
   nrmcheb(i) = log10(norm(yfit-y,inf));             
   concheb(i) = log10(cond(phi,2));   
end

% Plot approximation error per degree of approximation
figure
plot(n,nrmcheb,n,nrmunif)
xlabel('Degree of Approximating Polynomial'); ylabel('Log10 Error');
legend('Chebychev Nodes','Uniform Nodes','Location','N')
title('Log10 Polynomial Approximation Error for Runge''s Function');

% Plot interpolation matrix condition nnmbers per degree of approximation
figure
plot(n,concheb,n,conunif)
xlabel('Degree of Approximating Polynomial'); ylabel('Log10 Condition Number');
legend('Chebychev Polynomial Basis','Mononomial Basis','Location','N')
legend boxoff
title('Log10 Interpolation Matrix Condition Number');

% Plot Chebychev- and uniform node polynomial approximation errors
figure
hold on
plot(x,errcheb(:,5),x,errunif(:,5));
plot(x,0*x,'k--');
set(gca,'Xtick',[-1 0 1])
xlabel('x'); ylabel('Error');
title('Runge''s Function 11^{th}-Degree Polynomial Approximation Error.')
legend('Chebychev Nodes','Uniform Nodes','Location','N')

% Save Plots as EPS Files
printfigures(mfilename,4)


%% Function to be approximated
function y = runge(x)
y = 1./(1+25*x.^2);