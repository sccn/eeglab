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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
