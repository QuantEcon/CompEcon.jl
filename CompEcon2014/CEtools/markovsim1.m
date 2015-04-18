%% MARKOVSIM1
%
%  Markov chain simulator
%
%  Usage
%    j = markovsim(i,q)
%  Input
%    i   n.1  integers between 1 & m
%    q   m.m  probability transition matrix
%  Output
%    j = n.1  integers between 1 & m, randomly generated according to q

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu

function j = markovsim1(q,nr,np,i)

if nargin<4
  p = markov(q);
  i = discrand(nr,p);
  if np<2
    j = i;
    return
  end
end

m = size(q,1);
if np>1
  r = rand(nr,np-1);
  j = [i ones(nr,np-1)];
  for ip=2:np
    rr = r(:,ip-1);
    j(:,ip) = min(sum(rr(:,ones(1,m))>cumsum(q(j(:,ip-1),:),2),2)+1,m);
  end
else
  r = rand(nr,1);
  j = min(sum(r(:,ones(1,m))>cumsum(q(i,:),2),2)+1,m);
end