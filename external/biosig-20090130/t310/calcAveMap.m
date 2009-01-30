function r = calcAveMap(s, h, startend, trls, class, rmartif, montage, lowpass, cue, name)
% Displays the mean and standard deviation of each channel.
%
% This function calculates the mean and standard deviation of each channel and
% displays it in a topographical layout.
%
% Usage:
%   plotAveMap(s, h, startend, trls, class, rmartif, montage, lowpass, cue, name);
%
% Input parameters:
%   s        ... Input signal (as obtained by sload)
%   h        ... Header structure (as obtained by sload)
%   startend ... Start and end point within trial (s)
%
% Optional input parameters:
%   trls     ... Trials
%   class    ... Class to display
%   rmartif  ... Remove artifacts
%   montage  ... Electrode montage
%   lowpass  ... Lowpass filter cutoff frequency (Hz)
%   cue      ... Location of cue (s)
%   name     ... Name of the data set (or subject)
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
    error('Start and end point of one trial not specified.');
end;
if (nargin < 4)
    trls = [];
end;
if (nargin < 5)
    class = [];  % All classes
end;
if (nargin < 6)
    rmartif = true;
end;
if (nargin < 7)
    montage = [];
end;
if (nargin < 8)
    lowpass = [];
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

% Preprocess data
fs = h.SampleRate;

if (~isempty(lowpass))
    b = fir1(fs, lowpass/(fs/2));
    s = filtfilt(b, 1, s);
end;

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

% Topographic layout
if ~isempty(montage)
    if isnumeric(montage)
        plot_index = find(montage' == 1);
        n_rows = size(montage, 1);
        n_cols = size(montage, 2);
    elseif ischar(montage)
        [lap, plot_index, n_rows, n_cols] = getMontage(montage);
    end;
else
    plot_index = 1:size(s, 2);
    n_rows = ceil(sqrt(size(s, 2)));
    n_cols = n_rows;
end;

for k = 1:length(plot_index)
    temp = reshape(st(k,:)', triallen, trials)';

    r{k}.average = mean(temp);
    r{k}.stdev = std(temp);
    r{k}.cue = cue;
    r{k}.plot_index = plot_index;
    r{k}.n_rows = n_rows;
    r{k}.n_cols = n_cols;
    r{k}.startend = startend;
    r{k}.trials = trials;
    r{k}.class = class;
    r{k}.name = name;
    r{k}.fs = fs;
    r{k}.triallen = triallen;
    r{k}.chantype = h.CHANTYP(k);
    r{k}.datatype = 'Average map';
end;