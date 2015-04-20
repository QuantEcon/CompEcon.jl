function demlin02

%% DEMLIN02 Ill-conditioning of Vandermonde matrix

% Preliminary tasks
demosetup(mfilename)

% Compute approximation error and matrix condition number
warning off
n = (6:50)';
errv = zeros(length(n),1);
conv = zeros(length(n),1);
for i=1:length(n)
   v = vander(1:n(i)); 
   errv(i) = log10(norm(eye(n(i),n(i))-v\v)); 
   conv(i) = log10(cond(v));
end

% Smooth using quadratic function
X = [ones(size(n)) n];
b = X\errv
errv = X*b;
b = X\conv
conv= X*b;

% Plot matrix condition numbers
figure
plot(n,conv)
xlabel('n'); 
ylabel('Log_{10} Condition Number');
title('Vandermonde Matrix Condition Numbers');

% Plot approximation errors
figure
plot(n,errv)
xlabel('n'); 
ylabel('Log_{10} Error');
title('Approximation Error for I-V\\V')

% Save Plots as EPS Files
printfigures(mfilename,2)