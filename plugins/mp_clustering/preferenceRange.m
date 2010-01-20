%
% [pmin,pmax]=preferenceRange(s)
%
% Given a set of similarities, s, this function computes a lower
% bound, pmin, on the value for the preference where the optimal
% number of clusters (exemplars) changes from 1 to 2, and the
% exact value of the preference, pmax, where the optimal
% number of clusters changes from n-1 to n.
%
% For N data points, there may be as many as N^2-N pair-wise
% similarities (note that the similarity of data point i to k
% need not be equal to the similarity of data point k to i).
% These may be passed in an NxN matrix of similarities, s, where
% s(i,k) is the similarity of point i to point k. In fact, only
% a smaller number of relevant similarities need to be provided,
% in which case the others are assumed to be -Inf. M similarity
% values are known, can be passed in an Mx3 matrix s, where each
% row of s contains a pair of data point indices and a
% corresponding similarity value: s(j,3) is the similarity of
% data point s(j,1) to data point s(j,2).
%
% A single-cluster solution may not exist, in which case pmin is
% set to NaN.
%
% [pmin,pmax]=preferenceRange(s,METHOD) uses one of the methods
% below to compute pmin and pmax:
%
%   'exact'      Computes the exact values for pmin and pmax
%                (Warning: This can be quite slow)
%
%   'bound'      Computes the exact value for pmax, but estimates
%                pmin using a bound (default)
%
% Copyright (c) Brendan J. Frey and Delbert Dueck (2007). This
% software may be freely used and distributed for
% non-commercial purposes.

function [pmin,pmax]=preferenceRange(s,varargin);

% Check arguments
if nargin<1 error('Too few input arguments');
elseif nargin==1 method='bound';
elseif nargin==2 method=varargin{1};
else error('Too many input arguments');
end;

if length(size(s))~=2 error('s should be a 2D matrix');
elseif size(s,2)==3
    N=max(max(s(:,1)),max(s(:,2)));
    if min(min(s(:,1)),min(s(:,2)))<=0
        error('data point indices must be >= 1');
    end;
elseif size(s,1)==size(s,2)
    N=size(s,1);
else error('s must have 3 columns or be square'); end;

% Construct similarity matrix
if size(s,2)==3
    S=-Inf*ones(N,N); 
    for j=1:size(s,1) S(s(j,1),s(j,2))=s(j,3); end;
else S=s;
end;
for k=1:N S(k,k)=0; end;

% Find pmin
[dpsim1 k11]=max(sum(S,1));
if dpsim1==-Inf pmin=NaN;
elseif strcmp(method,'bound')
    for k=1:N S(k,k)=-Inf; end;
    m=max(S,[],2);
    tmp=sum(m);
    [yy ii]=min(m);
    tmp=tmp-yy-min(m([1:ii(1)-1,ii(1)+1:N]));
    pmin=dpsim1-tmp;
else
    dpsim2=-Inf;
    for j21=1:N-1
        for j22=j21+1:N
            tmp=sum(max(S(:,[j21,j22]),[],2));
            if tmp>dpsim2 dpsim2=tmp; k21=j21; k22=j22; end;
        end;
    end;
    pmin=dpsim1-dpsim2;
end;

% Find pmax
for k=1:N S(k,k)=-Inf; end;
pmax=max(S(:));
    