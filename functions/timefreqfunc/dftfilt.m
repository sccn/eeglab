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
end
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
 
