function r = calcHRr(s, h, triallen, trls, class, rmartif, ref, alpha, cue, name)
% Calculates and displays the instantaneous heart rate.
%
% This function calculates and displays the instantaneous heart rate relative to
% a reference interval. It also computes a confidence interval based on 
% bootstrap statistics.
%
% Usage:
%   [hr, cl, cu] = plotHRd(s, h, triallen, trls, class, rmartif, ref, alpha, cue, name);
%
% Input parameters:
%   s        ... Input signal (as obtained by sload)
%   h        ... Header structure (as obtained by sload)
%   triallen ... Length of one trial (s)
%
% Optional input parameters:
%   trls    ... Trials
%   class   ... Class
%   rmartif ... Remove artifacts
%   ref     ... Reference interval (in s)
%   alpha   ... Significance level
%   cue     ... Location of cue (s)
%   name    ... Name of the data set (or subject)
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
    error('Length of one trial not specified.');
end;
if (nargin < 4)
    trls = [];
end;
if (nargin < 5)
    class = [];  % All classes
end;
if (nargin < 6)
    rmartif = false;
end;
if (nargin < 7)
    ref = [];
end;
if (nargin < 8)
    alpha = [];
end;
if (nargin < 9)
    cue = [];
end;
if (nargin < 10)
    name = 'Subject';
end;

if (isempty(class))
    class = unique(h.Classlabel);
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

if (rmartif)
    [hrt, temp] = trigg(hr, h.TRIG(ismember(h.Classlabel, class) & ...
                        h.ArtifactSelection == 0), 0, triallen * h.SampleRate - 1);
else
    [hrt, temp] = trigg(hr, h.TRIG(ismember(h.Classlabel, class)), 0, ...
                        triallen * h.SampleRate - 1);
end;

% Select specified trials
if ~isempty(trls)
    hrt = reshape(hrt, temp);
    hrt = hrt(:, :, trls);
    hrt = reshape(hrt, temp(1), temp(2)*length(trls));
end;

triallen = triallen * fs;
trials = length(hrt) / triallen;

hrt = reshape(hrt', triallen, trials)';

if (isempty(ref))
    ref = mean(mean(hrt, 1), 2);  % Whole trial is reference
else
    ref = mean(mean(hrt(:,ref(1)*fs:ref(2)*fs), 1), 2);
end;

hrt = hrt/ref - 1;

if (~isempty(alpha))
    [hrt_boot, cl, cu] = bootts(hrt, 500, alpha);
    b = ones(1,5)/5;
    cl = filtfilt(b, 1, cl);
    cu = filtfilt(b, 1, cu);
end;

r{1}.cl = cl;
r{1}.cu = cu;
r{1}.hrt = hrt;
r{1}.triallen = triallen;
r{1}.trials = trials;
r{1}.alpha = alpha;
r{1}.cue = cue;
r{1}.name = name;
r{1}.class = class;
r{1}.fs = fs;
r{1}.ref = ref;
r{1}.datatype = 'Heart rate relative';