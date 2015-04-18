% CTBASEMAKE Basis matrices for continuous time collocation
% Computes the product of the coefficients and the basis matrices
% USAGE
%   B0=ctbasemake(rho,bases,0);
%   B1=ctbasemake(mu,bases,1);
% or
%   B2=ctbasemake(sigma,bases,2);
% In each case bases is created by a call to
%   bases=funbasex(basis,snodes,[0;1;2]);
% This creates a structured variable containing the primitive basis matrices.
% For d-dimensional state processes 
%   rho is N by 1
%   mu is N by d
% and
%   sigma is N by d by d or N by d^2
%
% If the first argument is omitted, a cell array of basis matrices is
% returned. This can be used in subsequent calls to ctbasemake. Thus
%    B1=ctbasemake(mu,bases,1);
% produces identical results as
%    Phi1=ctbasemake([],bases,1);
%    B1=ctbasemake(mu,Phi1,1);
% except that, with the latter, Phi1 may be reused.
% Using this option increases the memory demands but decreases computational
% time if B1 must be evaluated for multiple instances of mu.

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function B=ctbasemake(coeff,bases,order)
% bases are created by funbasex
if isa(bases,'struct')
  d=size(bases.order,2);
  [n,p]=size(coeff);
  switch order
  case 0
    B=bases.vals{1,d};
    for j=d-1:-1:1
      B=kron(B,bases.vals{1,j});
    end
    if ~isempty(coeff), 
      if p~=1, error('coeff must be nx1'); end
      B=spdiags(coeff,0,n,n)*B; 
    end
  case 1 
    index=eye(d)+1;
    if isempty(coeff)
      B=cell(1,d);
      for i=1:d
        B{i}=bases.vals{index(i,d),d};
        for j=d-1:-1:1
          B{i}=kron(B{i},bases.vals{index(i,j),j});
        end
      end
    else
      if p~=d, error('coeff must be nxd'); end
      B=[];
      for i=1:d
        if any(coeff(:,i)~=0)
          Bi=bases.vals{index(i,d),d};
          for j=d-1:-1:1
            Bi=kron(Bi,bases.vals{index(i,j),j});
          end
          if isempty(B), B=spdiags(coeff(:,i),0,n,n)*Bi;
          else, B=B+spdiags(coeff(:,i),0,n,n)*Bi;
          end
        end
      end
    end                    
  case 2
    u=ones(d,1); eyed=eye(d);
    order=kron(eyed,u)+kron(u,eyed);
    index=find(tril(ones(d,d))); % remove redundant evaluations
    index=order(index,:)+1;
    dd=d*(d+1)/2;
    if isempty(coeff)
      B=cell(1,dd);
      for i=1:dd
        B{i}=bases.vals{index(i,d),d};
        for j=d-1:-1:1
          B{i}=kron(B{i},bases.vals{index(i,j),j});
        end
      end
    else
      if prod(size(coeff))~=n*d*d, error('coeff must be nxdxd'); end
      coeff=rowvech(coeff,d);
      B=[];
      for i=1:dd
        if any(coeff(:,i)~=0)
          Bi=bases.vals{index(i,d),d};
          for j=d-1:-1:1
            Bi=kron(Bi,bases.vals{index(i,j),j});
          end
          if isempty(B), B=spdiags(coeff(:,i),0,n,n)*Bi;
          else, B=B+spdiags(coeff(:,i),0,n,n)*Bi;
          end
        end
      end
    end                    
  end
% bases are created by previous call to ctbasemake
else
  [n,p]=size(coeff);
  switch order
  case 0
    if p~=1, error('coeff must be nx1'); end
    B=spdiags(coeff,0,n,n)*bases;
  case 1
    d=size(bases,2); 
    if p~=d, error('coeff must be nxd'); end
    B=[];
    for i=1:d
      if isempty(B), B =   spdiags(coeff(:,i),0,n,n)*bases{i};
      else,          B = B+spdiags(coeff(:,i),0,n,n)*bases{i};
      end  
    end                    
  case 2
    d=round((1+sqrt(1+4*size(bases,2)))/2);
    if prod(size(coeff))/n~=d*d, error('coeff must be nxdxd'); end
    coeff=rowvech(coeff,d);
    B=[];
    for i=1:d*(d+1)/2
      if isempty(B), B =   spdiags(coeff(:,i),0,n,n)*bases{i};
      else,          B = B+spdiags(coeff(:,i),0,n,n)*bases{i};
      end
    end  
  end                  
end

if isempty(B), B=0; end



function Sigma=rowvech(sigma,d)
N=size(sigma,1);
Sigma=zeros(N,d*(d+1)/2);
sigma=reshape(sigma,N,d,d);
k=1;
for i=1:d
  Sigma(:,k)=sum(sigma(:,i,:).*sigma(:,i,:),3)/2;
  k=k+1;
  for j=i+1:d
    Sigma(:,k)=sum(sigma(:,i,:).*sigma(:,j,:),3);
    k=k+1;
  end
end
