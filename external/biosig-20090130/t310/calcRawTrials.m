function r = calcRawTrials(s, h, triallen, trls, class, rmartif, cue, name, plim)
% Prepares the data for displaying all trials of a specified class and channel.
%
% This function creates a struct that is used to display the raw data, arranged 
% in trials. Use the function plotC.m to plot this structure.
%
% Usage:
%   calcRawTrials(s, h, triallen, trls, class, rmartif, cue, name);
%
% Input parameters:
%   s        ... Input signal (as obtained by sload)
%   h        ... Header structure (as obtained by sload)
%   triallen ... Length of one trial (in s)
%
% Optional input parameters:
%   trls    ... Trials to display
%   class   ... Class
%   rmartif ... Remove artifacts
%   cue     ... Location of cue (in s)
%   name    ... Name of the data set (or subject)
%   plim    ... Scale y-axis to pmin for minimum and maximum value
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
    cue = [];
end;
if (nargin < 8)
    name = 'Subject';
end;
if (nargin < 9)
    plim = [];
end;

if (isempty(class))
    class = unique(h.Classlabel);
end;

% Trigger data
fs = h.SampleRate;

if (rmartif)
    [st, temp] = trigg(s, h.TRIG(ismember(h.Classlabel, class) & ...
                       h.ArtifactSelection == 0), 0, triallen * fs - 1);
else
    [st, temp] = trigg(s, h.TRIG(ismember(h.Classlabel, class)), 0, ...
                       triallen * fs - 1);
end;

% Select specified trials
if ~isempty(trls)
    st = reshape(st, temp);
    st = st(:, :, trls);
    st = reshape(st, temp(1), temp(2)*length(trls));
end;

triallen = triallen * fs;
trials = length(st) / triallen;

r{1}.trials = trials;
r{1}.triallen = triallen;
r{1}.st = st;
r{1}.cue = cue;
r{1}.class = class;
r{1}.fs = fs;
r{1}.plim = plim;
r{1}.name = name;
r{1}.datatype = 'Raw trials';