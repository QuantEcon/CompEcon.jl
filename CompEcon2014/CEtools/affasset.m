% AFFASSET
%
%  Solves affine asset pricing models
%
%  dX/dt = A'X + 0.5B'diag(C'X)C'X-g
%  dx/dt = a'X + 0.5b'diag(C'X)C'X-g0
%
% Usage
%   [X,x] = affasset(t,a,A,b,B,C,g,G,h,h0);
% Inputs
%   t              : m.1 vector of time values
%   a,A,b,B,C,g,g0 : parameters of n-state affine model
%   h,h0           : initial values h=X(t(1)) and h0=x(t(1)) 
% Ouutputs
%   X : m.n solution matrix for X(t)
%   x : m.1 solution vector for x(t)
%
% Uses: affode
%
% Note: CompEcon text refers to X and x as beta and beta0, respectively.

%  Copyright(c) 1997-2014
%    Paul L. Fackler  - paul_fackler@ncsu.edu
%    Mario J. Miranda - miranda.4@osu.edu

function [X,x] = affasset(t,a,A,b,B,C,g,g0,h,h0)

% Define parameters for ODE file (merges equations for X and x)
AA = [A.';a.'];
BB = [B.';b.']/2;
GG = [g(:);g0];

% Call solver
[t,X] = ode45('affode',t,[h(:);h0],[],AA,BB,C.',GG);

% Separate results
x = X(:,end);
X(:,end)=[];