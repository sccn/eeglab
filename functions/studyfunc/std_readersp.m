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
