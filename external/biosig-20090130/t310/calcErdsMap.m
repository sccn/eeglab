function r = calcErdsMap(s, h, startend, f, trls, class, rmartif, ref, submean, alpha, montage, cue, name)
% Calculates time-frequency (ERDS) maps.
%
% This function calculates time-frequency (ERDS) maps by using either FFT or
% wavelets to estimate the power in specific frequency bands. Maps can be
% calculated for more than one channel at once.
%
% Usage:
%   r = calcErdsMap(s, h, startend, f, trls, class, rmartif, ref, submean, alpha, montage, cue, name);
%
% Input parameters:
%   s        ... Input signal (as obtained by sload)
%   h        ... Header structure (as obtained by sload)
%   startend ... Start and end point within trial (s)
%   f        ... Frequency range (Hz)
%
% Optional input parameters:
%   trls    ... Trials to use (default: all)
%   class   ... Classes to use (default: all)
%   rmartif ... Remove artifacts (default: false)
%   ref     ... Reference interval (in s) (default: [] (whole trial))
%   submean ... Subtract mean signal first (default: yes)
%   alpha   ... Significance level (default: 0.05)
%   montage ... Electrode montage (default: [] (none))
%   cue     ... Location of cue (in s) (default: [] (none))
%   name    ... Name of the data set (or subject) (default: [] (none))
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
    error('Frequency range not specified.');
end;
if (nargin < 5)
    trls = [];  % All trials
end;
if (nargin < 6)
    class = [];  % All classes
end;
if (nargin < 7)
    rmartif = false;  % Don't remove artifacts
end;
if (nargin < 8)
    ref = [];  % Whole trial
end;
if (nargin < 9)
    submean = true;  % Subtract induced components
end;
if (nargin < 10)
    alpha = 0.05;  % Significance level of 5%
end;
if (nargin < 11)
    montage = [];  % No topographical layout
end;
if (nargin < 12)
    cue = [];  % No cue to draw
end;
if (nargin < 13)
    name = [];  % No subject name to draw
end;


if (isempty(class))
    class = unique(h.Classlabel);
end;

% Trigger data
fs = h.SampleRate;

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
cue = cue * 1000;  % The time needs to be in ms for erdsMap

% Subtract mean signal if requested
if (submean)
    for k = 1:size(s, 2)
        temp = reshape(st(k,:)', triallen, trials)';
        average = repmat(mean(temp), 1, trials);
        st(k, :) = st(k, :) - average;
    end;
end;

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
    n_cols = ceil(sqrt(size(s, 2)));
    if (size(s, 2) > 2)
        n_rows = n_cols;
    else
        n_rows = 1;
    end;
end;

for k = 1:size(s, 2)

    r{k} = erdsMap(st(k,:), triallen, [0 triallen/fs*1000], fs, ...
                   [round(f(1)) 0.75], 'maxfreq', round(f(2)), 'baseline', ...
                   ref.*1000, 'alpha', alpha, 'padratio', 8);
    r{k}.cue = cue;    
    r{k}.plot_index = plot_index;
    r{k}.n_rows = n_rows;
    r{k}.n_cols = n_cols;
    r{k}.triallen = triallen;
    r{k}.fs = fs;
    r{k}.f = f;
    r{k}.name = name;
    r{k}.trials = trials;
    r{k}.class = class;
    r{k}.triallen = triallen;
    r{k}.ref = ref;
    r{k}.alpha = alpha;
    r{k}.datatype = 'ERDS map';

end;