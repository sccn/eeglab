function bp = getBandpower(s, h, startend, f, trls, class, rmartif, submean)
% Calculates bandpower values in the specified frequency band.
%
% This function calculates bandpower values in the specified frequency band.
%
% Usage:
%   bp = getBandpower(s, h, triallen, f, trls, class, rmartif, submean);
% 
% Input parameters:
%   s        ... Input signal (as obtained by sload)
%   h        ... Header structure (as obtained by sload)
%   startend ... Start and end point within trial (s)
%   f        ... Frequency band (Hz)
%
% Optional input parameters:
%   trls    ... Trials
%   class   ... Class
%   rmartif ... Remove artifacts
%   submean ... Subtract mean signal
%
% Output parameters:
%   bp ... Bandpower

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
    error('Start and end point of one trial not specified.');
end;
if (nargin < 4)
    error('Frequency band not specified.');
end;
if (nargin < 5)
    trls = [];
end;
if (nargin < 6)
    class = [];  % All classes
end;
if (nargin < 7)
    rmartif = false;
end;
if (nargin < 8)
    submean = true;
end;

if (isempty(class))
    class = unique(h.Classlabel);
end;

% Preprocess data
fs = h.SampleRate;

[b, a] = butter(5, f/(fs/2));
s = filtfilt(b, a, s);

% Trigger data
if (rmartif)
    [st, temp] = trigg(s, h.TRIG(ismember(h.Classlabel, class) & ...
                       h.ArtifactSelection == 0), startend(1)*fs, ...
                       startend(2)*fs-1);
else
    [st, temp] = trigg(s, h.TRIG(ismember(h.Classlabel, class)), ...
                       startend(1)*fs, startend(2)*fs-1);
end;

% Select specified trials
if ~isempty(trls)
    st = reshape(st, temp);
    st = st(:, :, trls);
    st = reshape(st, temp(1), temp(2)*length(trls));
end;

triallen = (startend(2) - startend(1)) * fs;
trials = length(st) / triallen;

% Subtract mean signal if requested
if (submean)
    for k = 1:size(s, 2)
        temp = reshape(st(k,:)', triallen, trials)';
        average = repmat(mean(temp), 1, trials);
        st(k, :) = st(k, :) - average;
    end;
end;

% Calculate ERDS values
st = st.^2;

bp = zeros(size(s, 2), triallen);

for k = 1:size(s, 2)
    temp = reshape(st(k, :)', triallen, trials)';
    bp(k,:) = filtfilt(ones(1,round(fs/5))/(fs/5), 1, mean(temp));
end;