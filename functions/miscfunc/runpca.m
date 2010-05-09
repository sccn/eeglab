% runpca() -  perform principal component analysis (PCA) using singular value 
%             decomposition (SVD) using Matlab svd() or svds()
%                        >> inv(eigvec)*data = pc;
% Usage:
%    >> [pc,eigvec,sv] = runpca(data);
%    >> [pc,eigvec,sv] = runpca(data,num,norm)
%
% Inputs:
%   data   - input data matrix (rows are variables, columns observations)
%   num    - number of principal comps to return  {def|0|[] -> rows in data}
%   norm   - 1/0 = do/don't normalize the eigvec's to be equivariant 
%                                                {def|0 -> no normalization}
% Outputs:
%   pc     - the principal components, i.e.        >> inv(eigvec)*data = pc;
%   eigvec - the inverse weight matrix (=eigenvectors). >> data = eigvec*pc; 
%   sv     - the singular values (=eigenvalues)
%
% Author: Colin Humphries, CNL / Salk Institute, 1997
%
% See also: runica()

% Copyright (C) Colin Humphries, CNL / Salk Institute, Aug, 1997
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 01/31/00 renamed runpca() and improved usage message -sm
% 01-25-02 reformated help & license, added links -ad 

function [pc,M,S] = runpca(data,N,norm)

BIG_N = 50; % for efficiency, switch to sdvs() when BIG_N<=N or N==rows

if nargin < 1
  help runpca
  return
end

rows = size(data,1);

% remove the mean
for i = 1:rows
  data(i,:) = data(i,:) - mean(data(i,:));
end

if nargin < 3
  norm = 0;
elseif isempty(norm)
  norm = 0;
end

if nargin < 2
  N = 0;
end
if isempty(N)
  N = 0;
end


if N == 0  | N == rows
  N = rows;
  [U,S,V] = svd(data',0);   % performa SVD
  if norm == 0
    pc = U';
    M = (S*V')';
  else % norm
    pc = (U*S)';
    M = V;
  end
else
  if N > size(data,1)
    error('N must be <= the number of rows in data.')
  end
  %if N <= BIG_N | N == rows
  %[U,S,V] = svd(data',0);
  %else
  [U,S,V] = svds(data',N);
  %end
  if norm == 0
    pc = U';
    M = (S*V')';
  else % norm
    pc = (U*S)';
    M = V;
  end  
  %if N > BIG_N & N < rows
  %pc = pc(1:N,:);
  %M = M(:,1:N);
  %end
end
%S = diag(S(1:N,1:N));
