% spec() - power spectrum. This function replaces psd() function if the signal
%          processing toolbox is not present. It uses the timef() function.
%
% Usage:
%   >> [power freqs] = spec(X);
%   >> [power freqs] = spec(X, nfft, fs, win, overlap);
%
% Inputs:
%   X       - data
%   nfft    - zero padding to length nfft
%   fs      - sampling frequency
%   win     - window size
%   overlap - window overlap
%              
% Outputs:
%   power  - spectral estimate (amplitude not dB)
%   freqs  - frequency array
%
% Note: this function is just an approximation of the psd() (not pwelch) 
%       method. We strongly recommend to use the psd function if you have 
%       access to it.
%
% Known problems: 
%  1) normalization formula was determined manually by comparing
%     with the output of the psd function on about 100 examples.
%  2) In case only one time window is necessary, the overlapping factor 
%     will be increased so that at least 2 windows are presents (the 
%     timef function cannot use a single time window).
%  3) FOR FILTERED DATA, THE POWER OVER THE FILTERED REGION IS WRONG
%     (TOO HIGH)
%  4) the result of this function differs (in scale) from the pwelch 
%     function since the normalization is different for pwelch.
%
% Author: Arnaud Delorme, SCCN, Dec 2, 2003

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [power, freqs] = spec(X, nfft, fs, win, overlap);

if nargin < 1
    help spec;
    return;
end;

% default parameters
% ------------------
if nargin < 2
    nfft = 256;
else 
    nfft = pow2(nextpow2(nfft));
end;
nfft = min(length(X), nfft);
if nargin < 3
    fs = 2;
end;
if nargin < 4
    win = nfft;
else 
    win = pow2(nextpow2(win));
end;
if win > length(X)
    win = length(X);
end;
if log2(win) ~= round(log2(win))
    win = pow2(floor(log2(win)));
end;
if nargin < 5
    overlap = 0;
end;

% compute corresponding parameters for timef
% ------------------------------------------
padratio = pow2(nextpow2(nfft/win));
timesout = floor(length(X)/(win-overlap));
if timesout <= 1, timesout = 2; end;

[ersp itc mbase times freqs] = timef(X(:)', length(X), [0 length(X)]/fs*1000, fs, ...
                                        0, 'padratio', padratio, 'timesout', timesout, 'winsize', win, 'maxfreq', fs/2, ...
                                        'plotersp', 'off', 'plotitc', 'off', 'baseline', NaN, 'verbose', 'off');

ersp = 10.^(ersp/10); % back to amplitude
power = mean(ersp,2)*2.7/win; % this formula is a best approximation (I couldn't find the actual one)
                             % in practice the difference with psd is less than 0.1 dB
%power = 10*log10(power);
if nargout < 1
    hold on;
    h = plot(freqs,  10*log10(power));
    set(h, 'linewidth', 2);
end;
return;

figure;
stdv = std(ersp, [], 2);
h = plot(freqs, power+stdv, 'r:');
set(h, 'linewidth', 2);
h = plot(freqs, power-stdv, 'r:');
set(h, 'linewidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (log)');
