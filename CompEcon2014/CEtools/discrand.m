%% DISCRAND
%
%  Discrete random variable simulator
%
%  Usage
%    j = discrand(n,w)
%  Input
%    n = number of occurances simulated
%    w = m.1 vector of probabilities with which states 1,2,...,m occur
%  Output
%    j = n.1 vector of random occurances of states 1,2,...,m

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function j = discrand(n,w)
if n<1
  j = [];
else
  w = [0;cumsum(w(:))];
  m = length(w);
  m = m-length(find(w==w(end)));
  [~,j] = sort([w(1:m); rand(n,1)]);
  temp = find(j>m);
  i = j(temp)-m;
  j = temp-(1:n)';
  j(i) = j(:);
  j(j==0) = length(find(w==w(1)));
  j = max(j,1);
  j = min(j,m);
end