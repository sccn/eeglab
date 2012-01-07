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
%  'design'    - [integer] read files from a specific STUDY design. Default
%                is empty (use current design in STUDY.currentdesign).
%  'channels'  - [cell] list of channels to import {default: none}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if 
%                available). Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster
%                 {default:all}
%
% Spectrum specific inputs:
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'rmsubjmean' - ['on'|'off'] remove mean subject spectrum from every
%                 channel spectrum, making them easier to compare
%                 { default: 'off' }
% Output:
%  STUDY    - updated studyset structure
%  specdata - [cell array] spectral data (the cell array size is 
%             condition x groups)
%  freqs    - [float array] array of frequencies
%  setinds  - [cell array] datasets indices
%  cinds    - [cell array] channel or component indices
%
% Example:
%  std_precomp(STUDY, ALLEEG, { ALLEEG(1).chanlocs.labels }, 'spec', 'on');
%  [spec freqs] = std_readspec(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(1).labels });
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
    % old calling format
    % ------------------
    EEG = STUDY(ALLEEG);
    filename   = fullfile(EEG.filepath, EEG.filename(1:end-4));
    comporchan = varargin{1};
    options = {'measure', 'spec'};
    if length(varargin) > 1, options = { options{:} 'freqlimits', varargin{2} }; end;
    if comporchan(1) > 0
        [datavals tmp xvals] = std_readfile(filename, 'components',comporchan, options{:});
    else
       [datavals tmp xvals] = std_readfile(filename, 'channels', comporchan, options{:});
    end;
    STUDY    = datavals';
    specdata = xvals;
end;

[STUDY, specdata, allfreqs] = std_readerp(STUDY, ALLEEG, 'datatype', 'spec', varargin{:});
