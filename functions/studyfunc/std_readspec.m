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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [STUDY, specdata, allfreqs] = std_readspec(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readspec;
    return;
end

[STUDY, specdata, allfreqs] = std_readerp(STUDY, ALLEEG, 'datatype', 'spec', varargin{:});
