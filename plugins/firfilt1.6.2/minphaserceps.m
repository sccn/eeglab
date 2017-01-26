% rcepsminphase() - Convert FIR filter coefficient to minimum phase
%
% Usage:
%   >> b = minphaserceps(b);
%
% Inputs:
%   b - FIR filter coefficients
%
% Outputs:
%   bMinPhase - minimum phase FIR filter coefficients
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% References:
%   [1] Smith III, O. J. (2007). Introduction to Digital Filters with Audio
%       Applications. W3K Publishing. Retrieved Nov 11 2013, from
%       https://ccrma.stanford.edu/~jos/fp/Matlab_listing_mps_m.html
%   [2] Vetter, K. (2013, Nov 11). Long FIR filters with low latency.
%       Retrieved Nov 11 2013, from
%       http://www.katjaas.nl/minimumphase/minimumphase.html

function [bMinPhase] = minphaserceps(b)

% Line vector
b = b(:)';

n = length(b);
upsamplingFactor = 1e3; % Impulse response upsampling/zero padding to reduce time-aliasing
nFFT = 2^ceil(log2(n * upsamplingFactor)); % Power of 2
clipThresh = 1e-8; % -160 dB

% Spectrum
s = abs(fft(b, nFFT));
s(s < clipThresh) = clipThresh; % Clip spectrum to reduce time-aliasing

% Real cepstrum
c = real(ifft(log(s)));

% Fold
c = [c(1) [c(2:nFFT / 2) 0] + conj(c(nFFT:-1:nFFT / 2 + 1)) zeros(1, nFFT / 2 - 1)];

% Minimum phase
bMinPhase = real(ifft(exp(fft(c))));

% Remove zero-padding
bMinPhase = bMinPhase(1:n);
