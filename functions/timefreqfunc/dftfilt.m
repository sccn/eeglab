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
% Authors: Sigurd Enghoff, Arnaud Delorme & Scott Makeig, 
%          SCCN/INC/UCSD, La Jolla, 8/1/98

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

% 01-25-02 reformated help & license -ad 

% future developments
% -------------------
% input into dftfilt:
% - lowfreq and maxfreq (of interest)
% - lowcycle and maxcyle (ex: 3 cycles at low freq and 10 cycles at maxfreq)
% - the delta in frequency: ex 0.5 Hz
% The function should: compute the number of points (len) automatically
% Warning with FFT compatibility
% Still, something has to be done about the masking so that it would be comaptible

function b = dftfilt(len,maxfreq,cycle,oversmp,wavfact)

count = 1;
for index = 1:1/oversmp:maxfreq*len/cycle % scan frequencies
	w(:,count) = j * index * cycle * linspace(-pi+2*pi/len, pi-2*pi/len, len)'; % exp(-w) is a sinus curve
	count = count+1; % -2*pi/len ensures that we really scan from -pi to pi without redundance (-pi=+pi) 
end;
b = exp(-w);

%srate = 2*pi/len;						    % Angular increment.
%w = j * cycle * [0:srate:2*pi-srate/2]';	% Column.
%x = 1:1/oversmp:maxfreq*len/cycle;		    % Row.
%b = exp(-w*x);					            % Exponentiation of outer product.

for i = 1:size(b,2),
	m  = round(wavfact*len*(i-1)/(i+oversmp-1));	% Number of elements to discard.
	mu = round(m/2);				                % Number of upper elemnts.
	ml = m-round(m/2);				                % Number of lower elemnts.
	b(:,i) = b(:,i) .* [zeros(mu,1) ; hanning(len-m) ; zeros(ml,1)];
end

% syemtric hanning function
function w = hanning(n)
if ~rem(n,2)
   w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
   w = [w; w(end:-1:1)];
else
   w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
   w = [w; w(end-1:-1:1)];
end
 
