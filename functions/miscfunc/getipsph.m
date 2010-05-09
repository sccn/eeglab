% getipsph() - Compute "in place" (m by n) sphering or quasi-sphering matrix for an (n by t) 
%              input data matrix. Quasi-sphering reduces dimensionality of the data, while
%              maintaining approximately the "original" positions of the axes. That is, 
%              quasi-sphering "rotates back" as much as possible into the original channel 
%              axes, versus simple PCA reduction to the principal subspace (i.e., projecting 
%              the data onto a largest eigenvector basis).
% Usage:
%           >> S = getipsph(x,m)
%
% Input:
%           x  size [n,t] input n-channel data matrix
%
% Optional input:
%           m  (int <= n) Reduce output dimensions to [n,m]. If m is not present, or m==n, 
%              then the usual sphering matrix (the symmetric square root of the data 
%              covariance matrix) is returned. {default: [], return usual sphering matrix).
% Output:
%
%           S  size [m,n] (if m==n) sphering matrix, or (if m < n) quasi-sphering matrix 
% Example:
%           >> x = nrand(10,1000); % random (10,1000) matrix
%           >> S = getipsph(x,8);  % return quasi-sphering matrix reducing dimension to 8
%           >> d = S*x;            % d is the quasi-sphered 8-dimensional data
%           %
%           % If ICA decomposition is performed on d, as >> [w] = runica(d,'sphering','off');
%           % then the mixing matrix containing the 8 component maps in original channel 
%           % coordinates is >> A = pinv(W*S); % where A is size [n,m]
%

% Author: Jason Palmer, SCCN / INC / UCSD, 2008

% Copyright (C) Jason Palmer, SCCN / INC / UCSD , La Jolla 2008
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

function S = getipsph(x,m)

[n,N] = size(x);
if nargin < 2
    m = n;
end

if m>n
   help getipsph
   return
end

mn = mean(x,2);
for i = 1:n
    x(i,:) = x(i,:) - mn(i);
end
Sxx = x*x'/N; Sxx = (Sxx+Sxx')/2;
[U,D,V] = svd(Sxx);
ds = diag(D);

if m == n
    S = U * pinv(diag(sqrt(ds))) * U';
elseif m < n
    [sd,so] = sort(diag(Sxx),1,'descend');
    [uv,sv,vv] = svd(U(so(1:m),1:m));
    S = uv * vv' * pinv(diag(sqrt(ds(1:m)))) * U(:,1:m)';
else
    error('m must be less than or equal to n');
end

