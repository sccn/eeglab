function [EVec, Eload, EVal] = pca_calc(vecs);
%
%  PCA_CALC   Principal Component Analysis
%
%   [EVEC, ELOAD, EVAL] = PCA_CALC(X) takes a column-wise de-meaned
%       data matrix X of size n-by-m and calculates the n
%       Eigenvectors (EVEC of size n-by-n) of the data covariance
%       matrix, the factor loadings (ELOAD of size n-by-m) and the
%       corresponding Eigenvalues (EVAL of size n-by-1).
%
%
% Author: Rami K. Niazy, FMRIB Centre, University of Oxford
%
% Copyright (c) 2004 University of Oxford.
%

[m,n]=size(vecs);
[Us,S,EVec] = svd(vecs,0);

if m == 1
    S = S(1);
else
    S = diag(S);
end
Eload = Us .* repmat(S',m,1);
S = S ./ sqrt(m-1);   
if m <= n
    S(m:n,1) = 0;
    S(:,m:n) = 0;
end
EVal = S.^2;
