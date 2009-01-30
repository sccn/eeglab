function r = calcHRSpectrum(s, h)
% Calculates the spectrum of the heart rate.
%
% Calculates the heart rate spectrum of the signal s.
%
% Usage:
%   calcHRSpectrum(s, h);
%
% Input parameters:
%   s ... Input signal containing an ECG signal (as obtained by sload)
%   h ... Header structure (as obtained by sload)
%
% Output parameter:
%   r ... Structure containing the results

% Copyright by Robert Leeb, Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: robert.leeb@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Load data
hr_ges = [];


fs = h.SampleRate;

% Calculate heart rate
h_qrs = qrsdetect(s, fs, 2);
time = h_qrs.EVENT.POS;

% Correction of 3 samples
time = time - 3;

q = diff(time./fs);

bpm = 60./q;

% Plausibility check
for k = 2:length(bpm)
    if bpm(k) > (2*bpm(k-1))
        bpm(k) = bpm(k-1);
    elseif bpm(k) < (bpm(k-1)/2)
        bpm(k) = bpm(k-1);
    end;
end;

hr(1:time(1)) = bpm(1);
for k = 2:length(time)-1
    hr(time(k-1)+1:time(k)) = bpm(k);
end;

hr = hr(:);
hr(end+1:size(s,1)) = hr(end);
window = hanning(length(hr));
hr_ges = [hr_ges; hr.*window];


[p, f] = psd(hr_ges-mean(hr_ges), 200*fs, fs, hanning(100*fs), 50*fs, 'linear');
p = p / (50*fs);

r{1}.p = p;
r{1}.f = f;
r{1}.fs = fs;
r{1}.datatype = 'Heart rate spectrum';