function r = calcSpectrum(s, h, triallen, refint, actint, trls, class, rmartif, montage, f, name)
% Calculates the spectrum of each channel.
%
% This function calculates the spectrum of each channel in the input signal. Two
% conditions have to be specified, first the reference interval (to be plotted
% in blue), and second the activation interval (to be plotted in red).
%
% Usage:
%   r = calcSpectra(s, h, triallen, refint, actint, trls, class, rmartif, montage, f, name);
%
% Input parameters:
%   s        ... Input signal (as obtained by sload)
%   h        ... Header structure (as obtained by sload)
%   triallen ... Length of one trial (s)
%   refint   ... Reference interval (s)
%   actint   ... Activation interval (s)
%
% Optional input parameters:
%   trls    ... Trials
%   class   ... Class
%   rmartif ... Remove artifacts
%   montage ... Electrode montage
%   f       ... Upper frequency limit to display
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
    error('No reference interval specified.');
end;
if (nargin < 5)
    error('No activation interval specified.');
end;
if (nargin < 6)
    trls = [];
end;
if (nargin < 7)
    class = [];  % All classes
end;
if (nargin < 8)
    rmartif = false;
end;
if (nargin < 9)
    montage = [];
end;
if (nargin < 10)
    f = [];
end;
if (nargin < 11)
    name = 'Subject';
end;

if (isempty(class))
    class = unique(h.Classlabel);
end;
if (isempty(f))
    f = fs/2;
end;

% Trigger data
fs = h.SampleRate;

if (rmartif)
    [st_ref, temp_ref] = trigg(s, h.TRIG(ismember(h.Classlabel, class) & ...
                               h.ArtifactSelection == 0), refint(1) * fs, ...
                               refint(2) * fs - 1);
    [st_act, temp_act] = trigg(s, h.TRIG(ismember(h.Classlabel, class) & ...
                               h.ArtifactSelection == 0), actint(1) * fs, ...
                               actint(2) * fs - 1);
else
    [st_ref, temp_ref] = trigg(s, h.TRIG(ismember(h.Classlabel, class)), ...
                               refint(1) * fs, refint(2) * fs - 1);
    [st_act, temp_act] = trigg(s, h.TRIG(ismember(h.Classlabel, class)), ...
                               actint(1) * fs, actint(2) * fs - 1);
end;

% Select specified trials
if ~isempty(trls)
    st_ref = reshape(st_ref, temp_ref);
    st_ref = st_ref(:, :, trls);
    st_ref = reshape(st_ref, temp_ref(1), temp_ref(2)*length(trls));
    st_act = reshape(st_act, temp_act);
    st_act = st_act(:, :, trls);
    st_act = reshape(st_act, temp_act(1), temp_act(2)*length(trls));
end;

triallen = triallen * fs;
len_ref_int = (refint(2) - refint(1)) * fs;
len_act_int = (actint(2) - actint(1)) * fs;
trials = length(st_ref) / len_ref_int;

% Topographic layout
if ~isnan(montage)
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

r{1}.plot_index = plot_index;
r{1}.n_rows = n_rows;
r{1}.n_cols = n_cols;
r{1}.fs = fs;
r{1}.name = name;
r{1}.trials = trials;
r{1}.class = class;
r{1}.triallen = triallen;
r{1}.len_ref_int = len_ref_int;
r{1}.len_act_int = len_act_int;
r{1}.st_ref = st_ref;
r{1}.st_act = st_act;
r{1}.f = f;
r{1}.datatype = 'Spectrum';
