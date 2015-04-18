%% FUND
%
%  Evaluates functions and first 2 derivatives
%
%  Usage
%    [F,J,H] = fund(c,basis,x,HessOpt)
%  Input
%    c       : n.p basis coefficient matrix
%    basis   : approximation structure
%    x       : m.d evaluation points
%    HessOpt : 0/1/2 (default=0; see below)
%  Output
%    F : the function values
%    J : m.p.d Jacobian
%    H : Hessian
%  Note
%    J(i,j,k) is the derivative of the jth function with respect to the 
%    kth input variable evaluated at x(i,:).
%  Note
%    H(i,j,:) is the second derivative of the jth function with respect to 
%    input variables evaluated at x(i,:).  Form of H depends on HessOpt:
%        HessOpt = 0: vech form  - m.p.d(d+1)/2 array
%                = 1: vec form   - m.p.d^2 array
%                = 2: d x d form - m.p.d.d array

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [F,J,H] = fund(c,basis,x,HessOpt)

d = basis.d;
order = (0:nargout-1)';       % Only evaluate the bases that are needed
B = funbasex(basis,x,order);

% Evaluate the function
F = funeval(c,basis,B);
if nargout>1
  % Evaluate the Jacobian
  order = eye(d);
  J = funeval(c,basis,B,order);
  if nargout>2
    % Evaluate Hessian
    % Create the order matrix for the Hessian
    u = ones(d,1); 
    eyed = eye(d);
    order = kron(eyed,u)+kron(u,eyed);
    ind = tril(ones(d,d));
    order = order(ind,:);
    H = funeval(c,basis,B,order);
    % Expand the Hessian is vech form is not wanted
    if HessOpt>0
      H = H(:,:,vechinv(1:d*(d+1)/2));
      if HessOpt>1
        H = reshape(H,size(H,1),size(H,2),d,d);
      end
    end
  end
end