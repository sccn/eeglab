% dftfilt() - discrete Fourier filter
%
% Usage:
%   >> b = dftfilt(n,W,c,k,q)
%
% Inputs:
%   n - number of input samples
%   W - maximum angular freq. relative to n, 0 < W <= .5
%   c - cycles
%   k - oversampling
%   q - [0;1] 0->fft, 1->c cycles
%
% Authors: Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD, La Jolla, 8/1/98

% Copyright (C) 8/1/98 Sigurd Enghoff & Scott Makei, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license -ad 

function b = dftfilt(n,W,c,k,q)

f = 2*pi/n;						% Angular increment.
w = j * c * [0:f:2*pi-f/2]';	% Column.
x = 1:1/k:W*n/c;				% Row.
b = exp(-w*x);					% Exponentiation of outer product.

for i = 1:size(b,2),
	m  = round(q*n*(i-1)/(i+k-1));	% Number of elements to discard.
	mu = round(m/2);				% Number of upper elemnts.
	ml = m-round(m/2);				% Number of lower elemnts.
	b(:,i) = b(:,i) .* [zeros(mu,1) ; hanning(n-m) ; zeros(ml,1)];
%	b(:,i) = b(:,i) .* [zeros(mu,1) ; ones(n-m,1) ; zeros(ml,1)];
end
