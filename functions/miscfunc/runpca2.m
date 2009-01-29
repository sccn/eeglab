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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.6  2007/05/12 02:36:14  jason
% added code to remove mean, which pca should do.
%
%
%
%
% ls
% ls
% pwd
%
% Revision 1.5  2004/05/18 18:02:15  arno
% commenting 2 not usefull lines
%
% Revision 1.4  2004/05/05 16:30:57  arno
% svd -> svds
%
% Revision 1.3  2004/05/05 15:22:20  arno
% remove svds
%
% Revision 1.2  2004/05/05 01:37:17  arno
% lowercase n -> uppercase N
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 01/31/00 renamed runpca() and improved usage message -sm
% 01-25-02 reformated help & license, added links -ad 

function [pc,sph]= runpca2(data,N)

if nargin < 1
  help runpca
  return
end

[chans,frames] = size(data);

% remove the mean
for i = 1:chans
  data(i,:) = data(i,:) - mean(data(i,:));
end

% do svd
[U,S,V] = svd(data*data'/frames); % U and V should be the same since data*data' is symmetric

sv = diag(S);
pc = U(:,1:N) * diag(sv(1:N)) * V(:,1:N)';
sph = U(:,1:N) * pinv(diag(sqrt(sv(1:N)))) * V(:,1:N)';
