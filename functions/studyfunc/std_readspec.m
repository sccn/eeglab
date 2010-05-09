% std_readspec() - load spectrum measures for data channels or 
%                  for all components of a specified cluster.
%                  Called by plotting functions
%                  std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> [STUDY, specdata, allfreqs, setinds, cinds] = ...
%                   std_readspec(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'channels'  - [cell] list of channels to import {default: all}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'rmsubjmean' - ['on'|'off'] remove mean subject spectrum from every
%                 channel spectrum, making them easier to compare
%                 { default: 'off' }
%  'subject'    - [string] select a specific subject {default:all}
%  'component'  - [integer] select a specific component in a cluster
%                 {default:all}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if available)
%
% Output:
%  STUDY    - updated studyset structure
%  specdata - [cell array] spectral data (the cell array size is 
%             condition x groups)
%  freqs    - [float array] array of frequencies
%  setinds  - [cell array] datasets indices
%  cinds    - [cell array] channel or component indices
%
% Author: Arnaud Delorme, CERCO, 2006-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [STUDY, specdata, allfreqs] = std_readspec(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readspec;
    return;
end
if ~isstruct(ALLEEG) % old calling format
    [STUDY, specdata, allfreqs] = std_readspecsub(STUDY, ALLEEG, varargin{:});
    return;
end;

[STUDY, specdata, allfreqs] = std_readerp(STUDY, ALLEEG, 'datatype', 'spec', varargin{:});
return;

% std_readspecsub() - returns the stored mean power spectrum for an ICA component 
%                  in a specified dataset.  The spectrum is assumed to have been 
%                  saved in a Matlab file, "[dataset_name].icaspec", in the same
%                  directory as the dataset file. If this file doesn't exist,
%                  use std_spec() to create it or a pre-clustering function
%                  (pop_preclust() or std_preclust()) that calls it. 
% Usage:    
%  >> [spec, freqs] = std_readspecsub(ALLEEG, setindx, component, freqrange, rmsubjmean);  
%
% Inputs:
%   ALLEEG     - a vector of dataset EEG structures (may also be one dataset). 
%                Must contain the dataset of interest (the 'setindx' below).
%   setindx    - [integer] an index of an EEG dataset in the ALLEEG
%                structure for which to read a component spectrum.
%   component  - [integer] index of the component in the selected EEG dataset 
%                for which to return the spectrum
%   freqrange  - [min max in Hz] frequency range to return
%   rmsubjmean - [0|1] remove subject mean spectrum (0 is no and is the default)
%
% Outputs:
%   spec      - the log-power spectrum of the requested ICA component in the
%               specified dataset (in dB)
%   freqs     - vector of spectral frequencies (in Hz)
%
%  See also  std_spec(), pop_preclust(), std_preclust()
%
% Authors:  Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

function [X, f, singletrialdatapresent] = std_readspecsub(setinfo, chancomp, freqrange, rmsubjmean, singletrial)

if nargin < 3
    freqrange = [];
end;
if nargin < 4
    rmsubjmean = 0;
end;
if nargin < 5
    singletrial = 0;
end;

X = [];
f = [];
chanlab  = {};
if iscell(chancomp)
    if isfield(setinfo,'filebase')
        locs = load('-mat', [ setinfo(1).filebase '.datspec'], 'labels');
        if ~isempty(locs)
            chan.chanlocs = struct('labels', locs.labels);
            chancomp = -std_chaninds(chan, chancomp);
        else
            warning('Recomputing data file, old version')
            return;
        end;
    elseif isfield(setinfo,'chanlocs')
        chans = setinfo;
        chancomp = -std_chaninds(chans, chancomp);
    end;
end;
if chancomp(1) < 0
    chanorcomp = 'chan';
else
    chanorcomp = 'comp';
end;
if length(chancomp) < length(setinfo)
    chancomp(1:length(setinfo)) = chancomp(1);
end;

singletrialdatapresent = 1;
for k = 1:length(setinfo)

    % convert chancomponents or channel indices
    % -------------------------------------
    if iscell(chancomp)
        error('Cannot process cell array');
    elseif chancomp(1) < 0
        inds   = -chancomp;
    else
        inds   = chancomp;
    end;

    if strcmpi(chanorcomp, 'chan')
         filename = [ setinfo(k).filebase '.datspec'];
    else filename = [ setinfo(k).filebase '.icaspec'];
    end;

    %try,
        warning('off', 'MATLAB:load:variableNotFound');
        if rmsubjmean == 0
             erpstruct = load( '-mat', filename, [ chanorcomp int2str(inds(k)) ], 'freqs' );
        else erpstruct = load( '-mat', filename, [ chanorcomp int2str(inds(k)) ], 'freqs', 'average_spec' );
        end;
        warning('on', 'MATLAB:load:variableNotFound');
    %catch
    %    error( [ 'Cannot read file ''' filename '''' ]);
    %end;

    tmpdat    = getfield(erpstruct, [ chanorcomp int2str(inds(k)) ]);
    if singletrial == 0,
        if size(tmpdat,2) > 1 && size(tmpdat,1) > 1, tmpdat = mean(tmpdat,2); end;
    else
        if size(tmpdat,1) == 1 || size(tmpdat,2) == 1
            singletrialdatapresent = 0;
        end;
    end;
    if rmsubjmean, 
        if isfield(erpstruct, 'average_spec')
            if singletrial == 0
                 tmpdat = tmpdat - erpstruct.average_spec; 
            else tmpdat = tmpdat - repmat(erpstruct.average_spec', [1 size(tmpdat,2)]); 
            end;
        end;
    end;
    if singletrial
        if k == 1
            if size(tmpdat,1) == 1, tmpdat = tmpdat'; end;
            X = zeros([ 1 size(tmpdat) ]);
            X(1,:,:) = tmpdat;
        else
            X(1,:,end+1:end+size(tmpdat,2)) = tmpdat;
        end;
    else
        if k == 1
            if size(tmpdat,1) == 1, tmpdat = tmpdat'; end;
            X = zeros([ length(chancomp) size(tmpdat) ]);
        end;
        X(k,:,:) = tmpdat;
    end;
    f = getfield(erpstruct, 'freqs');
end;

% select frequency range of interest
% ----------------------------------
if ~isempty(freqrange)
    maxind = max(find(f <= freqrange(end)));
    minind = min(find(f >= freqrange(1)));
else
    %if not, use whole spectrum
    maxind = length(f);
    minind = 1;
end

f = f(minind:maxind);
X = X(:,minind:maxind,:);
return;
