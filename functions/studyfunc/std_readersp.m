% std_readersp() - load ERSP measures for data channels or  for all 
%                  components of a specified cluster. This function is 
%                  also being used to read ITC and ERPimage data.
% Usage:
%   >> [STUDY, erspdata, times, freqs, erspbase] = ...
%                   std_readersp(STUDY, ALLEEG, varargin);
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
%  'forceread' - ['on'|'off'] Force rereading data from disk.
%                Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster.
%                This is the index of the component in the cluster not the
%                component number {default:all}
%  'datatype'  - {'ersp'|'itc'|'erpim'} This function is used to read all 
%                2-D STUDY matrices stored on disk (not only ERSP). It may
%                read ERSP ('ersp' option), ITC ('itc' option) or ERPimage
%                data ('erpim' option).
%
% ERSP specific options:
%  'timerange' - [min max] time range {default: whole measure range}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'subbaseline' - ['on'|'off'] subtract the ERSP baseline for paired
%                conditions. The conditions for which baseline is removed
%                are indicated on the command line. See help
%                std_studydesign for more information about paired and
%                unpaired variables.
%
% ERPimage specific option:
%   This function is used to read all 2-D STUDY matrices stored on disk
%   (this includes ERPimages). It therefore takes as input specific 
%   ERPimage options. Note that the 'singletrials' optional input is
%   irrelevant for ERPimages (which are always stored as single trials).
%   'concatenate' - ['on'|'off'] read concatenated ERPimage data ('on') or 
%                   stacked ERPimage data. See help std_erpimage for more 
%                   information.
%   'timerange'   - [min max] time range {default: whole measure range}
%   'trialrange'  - [min max] read only a specific range of the ERPimage
%                   output trials {default: whole measure range}
% Output:
%  STUDY    - updated studyset structure
%  erspdata - [cell array] ERSP data (the cell array size is 
%             condition x groups). This may also be ITC data or ERPimage
%             data (see above).
%  times    - [float array] array of time points
%  freqs    - [float array] array of frequencies. For ERPimage this
%             contains trial indices.
%  erspbase - [cell array] baseline values
%  events   - [cell array] events (ERPimage only).
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

function [STUDY, erspdata, alltimes, allfreqs, erspbase, events, unitPower] = std_readersp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readersp;
    return;
end
events = {};
unitPower = 'dB';
erspbase = [];

disp('This function is obsolete and only partially backward compatible use std_readdata instead');
[STUDY, erspdata, alltimes, allfreqs, events] = std_readdata(STUDY, ALLEEG, 'datatype', 'ersp', varargin{:});
