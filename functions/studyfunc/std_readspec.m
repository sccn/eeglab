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

% $Log: std_readspec.m,v $
% Revision 1.20  2010/02/24 15:20:22  claire
% typo for error message
%
% Revision 1.19  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%
% Revision 1.18  2007/02/05 16:17:13  arno
% fix crash for old study subtracting average spectrum
%
% Revision 1.17  2006/11/23 01:12:50  arno
% implement mean spectrum subtraction
%
% Revision 1.16  2006/11/10 00:08:43  arno
% reprogram for channel
%
% Revision 1.15  2006/10/04 23:39:49  toby
% Bug fix courtesy Bas de Kruif
%
% Revision 1.14  2006/03/28 15:38:13  scott
% help msg
%
% Revision 1.13  2006/03/14 02:32:32  scott
% help msg
%
% Revision 1.12  2006/03/11 07:30:01  arno
% freqrange input
%
% Revision 1.11  2006/03/11 07:25:37  arno
% header
%
% Revision 1.10  2006/03/10 16:33:37  arno
% selecting frequency range for reading
%
% Revision 1.9  2006/03/10 00:37:45  arno
% error msg
%
% Revision 1.8  2006/03/09 18:10:38  arno
% *** empty log message ***
%
% Revision 1.7  2006/03/09 18:10:18  arno
% do not use etc field any more
%
% Revision 1.6  2006/03/09 00:42:09  arno
% fix reading file
%
% Revision 1.5  2006/03/09 00:37:31  arno
% now writing matlab fileend
%
% Revision 1.4  2006/03/09 00:03:57  arno
% read spectrum form matlab file
%
% Revision 1.3  2006/03/08 21:06:37  arno
% rename func
%
% Revision 1.2  2006/03/07 22:21:12  arno
% use fullfile
%

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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: std_readspec.m,v $
% Revision 1.20  2010/02/24 15:20:22  claire
% typo for error message
%
% Revision 1.19  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%

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
