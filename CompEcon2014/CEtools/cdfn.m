% CDFN Computes the CDF of the standard normal distribution
% USAGE
%   p=cdfn(x)
% A C-mex file version is available
function p=cdfn(x)
p = 0.5 * erfc(-0.7071067811865475*x);   % x/-sqrt(2)
