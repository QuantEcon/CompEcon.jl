%% MARKOVSIM
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

function j = markovsim(i,q)
m = size(q,1);
n = length(i);
r = rand(n,1);
j = min(sum(r(:,ones(1,m))>cumsum(q(i,:),2),2)+1,m);