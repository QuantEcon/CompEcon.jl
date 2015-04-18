% CHEBBASEX A utility used by CHEBBASE
% Coded as a C-mex file for additional speed

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
function bas=chebbasex(n,a,b,x)
z = (2/(b-a))*(x-(a+b)/2);
m=size(z,1);
bas=zeros(m,n);
bas(:,1)=ones(m,1);
bas(:,2)=z;
z=2*z;
for i=3:n
  bas(:,i)=z.*bas(:,i-1)-bas(:,i-2);
end