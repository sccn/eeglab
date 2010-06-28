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

function [pc,A,sv]= runpca2(data,K,dosym)
[chans,frames] = size(data);

if nargin < 3 || K < chans
    dosym = 1;
end
if nargin < 2
    K = chans;
end
if nargin < 1
  help runpca
  return
end

% remove the mean
for i = 1:chans
  data(i,:) = data(i,:) - mean(data(i,:));
end

% do svd
[U,S,V] = svd(data*data'/frames); % U and V should be the same since data*data' is symmetric

sv = sqrt(diag(S));
if dosym == 1
    pc = pinv(diag(sv(1:K))) * V(:,1:K)' * data;
    A = U(:,1:K) * diag(sv(1:K));
else
    pc = U * pinv(diag(sv)) * V' * data;
    A = U * diag(sv) * V';
end
