% QZORDERED QZ decomposition ordered by the absolute value of the generalized eigenvalues
% See QZ

% Based on code by Pedro Oviedo & Chris Sims

function [S,T,Q,Z]=qzordered(A,B);
n = size(A,1);		
[S,T,Q,Z]=qz(A,B);
i = 1;
while i<n;
   if abs(T(i,i)*S(i+1,i+1))>abs(S(i,i)*T(i+1,i+1));    
      [S,T,Q,Z] = qzswitch(i,S,T,Q,Z);
      if (i>1), i = i-1; 
      else,     i = i+1; 
      end
   else
     i=i+1;
   end
end



function [S,T,Q,Z] = qzswitch(i,S,T,Q,Z)

a = S(i,i); d = T(i,i); b = S(i,i+1); e = T(i,i+1);
c = S(i+1,i+1); f = T(i+1,i+1); 
wz = [c*e-f*b, (c*d-f*a)'];
n = sqrt(wz*wz');

if n == 0
   return
else
   xy = [(b*d-e*a)', (c*d-f*a)'];
   m = sqrt(xy*xy');
   wz = wz/n;
   xy = xy/m;
   wz = [wz; -wz(2)', wz(1)'];
   xy = [xy; -xy(2)', xy(1)'];
   S(i:i+1,:) = xy*S(i:i+1,:);
   T(i:i+1,:) = xy*T(i:i+1,:);
   S(:,i:i+1) = S(:,i:i+1)*wz;
   T(:,i:i+1) = T(:,i:i+1)*wz;
   Z(:,i:i+1) = Z(:,i:i+1)*wz;
   Q(i:i+1,:) = xy*Q(i:i+1,:);
end
