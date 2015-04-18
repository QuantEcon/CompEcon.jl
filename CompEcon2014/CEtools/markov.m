% MARKOV Analyzes Markov transition probability matrices
% USAGE
%   [q,f]=markov(p)
% INPUT
%   P :  nxn Markov transition probability matrix
%        (non-negative matrix with unit row sums)
% OUTPUTS
%    Q : nxk matrix of invariant distributions 
%          for each recurrence class
%        Q(i,k)=long-run frequency of visits to state i
%          given the system enters recurrence class k
%    F : nxn matrix of accessibility probabilities
%        F(i,j)= Prob[state j will be reached from state i]
%
% Useful in analyzing long-run characteristics of the optimized
% transition probability matrix resulting from the solution of discrete
% dynamic programming problems, e.g., the matrix PSTAR returned by
% DDPSOLVE.
%
% Reference: "Introduction to Stochastic Processes" by Ethan Cinlar
%   Prentice-Hall, 1975. Chapters 5 and 6.
%
% See also: DDPSOLVE.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [q,f]=markov(p);

n=size(p,1);

% Error Checking to ensure P is a valid stochastic matrix
if size(p,2)~=n
  error('Transition matrix is not square');
end
if any(any(p<0))
  error('Transition matrix contains negative elements');
end
if any(abs(p*ones(n,1)-1)>1e-12)
  error('Rows of transition matrix do not sum to 1');
end


% Determine accessibility from i to j
f=zeros(n,n);
for j=1:n
  dr=1;
  r=spones(p(:,j));            % a vector of ones where p(i,j)~=0
  while any(dr)
    dr=r;
    r=spones(p*r+r);
    dr=r-dr;
  end
  f(:,j)=r;
end


% Determine membership in recurrence classes
% Programming note:
%  f(i,:)=1 for states accessible from i
%  f(:,i)'=1 for states from which i is accessible
%  f(:,i)'.*f(i,:)=1 for states communicating with i (two-way accessibility)
%  If the set of communicating states is the same as the set of accessible 
%    states, it forms a recurrence class.
ind=zeros(n,n);
numrec=0;                          % number of recurrence classes
for i=1:n
  if all(ind(i,:)==0)
    j=f(i,:);               % states accessible from i
    if all((f(:,i)'.*j)==j) % are all accessible states communicating states?
      j=find(j);            % members in class with state i
      k=length(j);          % # of members
      if k>0 
        numrec=numrec+1; 
        ind(j,numrec)=ones(k,1); 
      end
    end
  end
end
ind=ind(:,1:numrec);        % ind(i,j)=1 if state i is in class j

% Determine recurrence class invariant probabilities
q=zeros(n,numrec);
for j=1:numrec
  k=find(ind(:,j));             % members in class j
  nk=length(k);                 % # of members
  % solve Pq=q s.t. 1'q=1
  q(k,j)=[ones(1,nk);(speye(nk)-p(k,k)')]\[1;zeros(nk,1)];
end

% Analyze transients if second output desired
if nargout>1
if numrec>1 trans=find(sum(ind')==0);
else trans=find(ind==0);
end
numt=length(trans);          % number of transient states
% Determine transient absorption and reachability probabilities
if numt>0                    
  qq=p(trans,trans);
  b=zeros(numt,n);
  for j=1:numrec
    k=find(ind(:,j));        % members of class j
    nk=length(k);            % # of members
    if nk==1                 % 1-step prob: transient states to class j
      b(:,k)=p(trans,k);
    else
      b(:,k)=sum(p(trans,k)')'*ones(1,nk);
    end
  end;
  qq=inv(eye(numt)-qq);          
  f(trans,:)=qq*b;              % absorption probabilities
  d=diag(qq)';
  qq=qq./d(ones(numt,1),:);
  f(trans,trans)=qq-diag(1./d); % transient reachability probabilities
end
end
