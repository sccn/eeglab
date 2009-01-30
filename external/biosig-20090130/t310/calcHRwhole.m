function r = calcHRwhole(s, h, startend)
% Calculates the instantaneous heart rate.
%
% This function calculates the instantaneous heart rate. It is also possible to 
% mark the occurrences of cues by vertical lines.
%
% Usage:
%   hr = plotHR(s, h, startend);
%
% Input parameters:
%   s ... Input signal (as obtained by sload)
%   h ... Header structure (as obtained by sload)
%
% Optional input parameters:
%   startend ... Start and end of the signal to be plotted (s)
%
% Output parameter:
%   r ... Structure containing the results

% Copyright by Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: clemens.brunner@tugraz.at

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

if (nargin < 1)
    error('No input signal specified.');
end;
if (nargin < 2)
    error('No header structure specified.');
end;
if (nargin < 3)
    startend = [];  % Plot whole signal
end;

% Calculate heart rate
fs = h.SampleRate;

if (size(s, 2) ~= 1)
    chanecg = find(h.CHANTYP == 'C');
    chanecg = chanecg(1);
else
    chanecg = 1;
end;

h_qrs = qrsdetect(s(:, chanecg), fs, 2);
time = h_qrs.EVENT.POS;

% Correction of 3 samples
time = time - 3;

q = diff(time./fs);

bpm = 60./q;

% Plausibility check
for i = 2:length(bpm)
    if bpm(i) > (2*bpm(i-1))
        bpm(i) = bpm(i-1);
    elseif bpm(i) < (bpm(i-1)/2)
        bpm(i) = bpm(i-1);
    end;
end;

hr(1:time(1)) = bpm(1);
for j = 2:length(time)-1
    hr(time(j-1)+1:time(j)) = bpm(j);
end;

hr = hr(:);
hr(end+1:size(s,1)) = hr(end);

r{1}.hr = hr;
r{1}.fs = fs;
r{1}.startend = startend;
r{1}.h = h;
r{1}.datatype = 'Heart rate continuous';