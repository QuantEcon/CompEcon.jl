% RSRES residual function for regime switching solver
% Used by RSSOLVE

function [e,c]=rsres(y,x,func,params,xindex,xchoice,m,n,...
              basis,z,ncond,Phi0,Phi1,Phi2,phi0,phi1,phi2)

  nm=sum(n);
  x(xchoice(:))=y; 
  for i=1:m
    j=xindex(:,1)==i | xindex(:,2)==i;
    a(i)=min(x(j));
    b(i)=max(x(j));
  end
  r=b-a;
  
  % loop over the regimes
  H=zeros(nm,nm);
  h=zeros(nm,1);
  rows=0;
  cols=0;
  for i=1:m
    S=z{i}*r(i)+a(i);
    nz=n(i)-ncond(i);
    rows=rows(end)+1:rows(end)+nz;
    cols=cols(end)+1:cols(end)+n(i);
    rho=feval(func,'rho',S,i,params{:});
    mu=feval(func,'g',S,i,params{:});
    sigma=feval(func,'sigma',S,i,params{:}); 
    if ~isempty(sigma)
       H(rows,cols)=spdiags(rho,0,nz,nz)*Phi0{i}...
          - spdiags(mu/r(i),0,nz,nz)*Phi1{i}...
          - spdiags(sigma.*sigma/(2*r(i)^2),0,nz,nz)*Phi2{i};
    else
       H(rows,cols)=spdiags(rho,0,nz,nz)*Phi0{i}...
          -spdiags(mu/r(i),0,nz,nz)*Phi1{i};
    end
    h(rows)=feval(func,'f',S,i,params{:});
  end
  
  % Loop over the switch/boundary points
  R=feval(func,'reward',x,[],params{:});
  G=[]; g=[];
  row=rows(end);
  for k=1:size(x,1)
    i=xindex(k,1);
    zi=(x(k)-a(i))/r(i);
    icols=sum(n(1:i-1))+(1:n(i));
    j=xindex(k,2);
    if j~=0
      zj=(x(k)-a(j))/r(j); 
      jcols=sum(n(1:j-1))+(1:n(j));
    end
    for order=0:2
      if xindex(k,3+order)==1 
        row=row+1;  
        H(row,icols)=funbase(basis{i},zi,order)/r(i)^order;
        if j~=0
          H(row,jcols)=-funbase(basis{j},zj,order)/r(j)^order;
        end
        h(row)=R(k,order+1);
      end
      if xindex(k,3+order)==2
        temp=zeros(1,nm);
        temp(1,icols)=funbase(basis{i},zi,order)/r(i)^order;
        if j~=0
          temp(1,jcols)=-funbase(basis{j},zj,order)/r(j)^order;
        end
        G=[G;temp];
        g=[g;R(k,order+1)];
      end
    end
  end

  c=H\h;
  e=G*c-g;

return