% std_readcustom() - Read custom data structure for file save on disk.
%
% Usage:    
% >> data = std_readcustom(STUDY, ALLEEG, fileext, 'key', 'val', ...);
%
% Required inputs:
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   fileext      - [string] file extension (without '.')
%
% Optional inputs:
%  'design'    - [integer] use specific study index design to compute measure.
%                Default is to use the default design.
%  'datafield' - [string or cell] extract only specific variables from the data
%                files. By default, all fields are loaded. Use '*' to match
%                patterns. If more than 1 variable is selected, data is
%                placed in a structure named data.
%  'eegfield'  - [string] copy data to a field of an EEG structure and return
%                EEG structure. Default is to return the data itself.
%  'eegrmdata' - ['on'|'off'] when option above is used, remove data from 
%                EEG structures before returning them. Default is 'on'.
% 
% Outputs:
%   data - cell array containing data organized according to the selected
%          design.
%
% Example:
%   % assuming ERP have been computed for the currently selected design
%   data = std_readcustom(STUDY, ALLEEG, 'daterp', 'datafield', 'chan1');
%   data = cellfun(@(x)x', siftdata, 'uniformoutput', false); % transpose data
%   std_plotcurve([1:size(data{1})], data); % plot data
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2013-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2013, arno@sccn.ucsd.edu
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

function [ returndata ] = std_readcustom(STUDY, ALLEEG, fileext, varargin)
    
    if nargin < 2
        help std_sift;
        return;
    end
    
    [g arguments] = finputcheck(varargin, { 'design'     'integer'           [] STUDY.currentdesign;
                                            'datafield'  { 'string' 'cell' } [] {};
                                            'eegfield'   'string'            [] '';
                                            'eegrmdata'  'string'            { 'on' 'off' } 'on' }, 'std_sift', 'mode', 'ignore');
    if ischar(g), error(g); end
    if ~iscell(g.datafield), g.datafield = { g.datafield }; end
    
    % Scan design and save data
    % -------------------------
    nc = max(length(STUDY.design(g.design).variable(1).value),1);
    ng = max(length(STUDY.design(g.design).variable(2).value),1);
    
    for cInd = 1:nc
        for gInd = 1:ng
            
            % find the right cell in the design
            cellInds = [];
            for index = 1:length(STUDY.design(g.design).cell)
                condind = std_indvarmatch( STUDY.design(g.design).cell(index).value{1}, STUDY.design(g.design).variable(1).value);
                grpind  = std_indvarmatch( STUDY.design(g.design).cell(index).value{2}, STUDY.design(g.design).variable(2).value);
                if isempty(STUDY.design(g.design).variable(1).value), condind = 1; end
                if isempty(STUDY.design(g.design).variable(2).value), grpind  = 1; end
                if cInd == condind && gInd == grpind
                    cellInds = [ cellInds index ];
                end
            end
            
            desset = STUDY.design(g.design).cell(cellInds);
            clear EEGTMP data;
            for iDes = 1:length(desset)
                
                % load data on disk
                tmpData = load('-mat', [ STUDY.design(g.design).cell(cellInds(iDes)).filebase '.' fileext ], g.datafield{:});
                
                % put data in EEG structure if necessary
                if ~isempty(g.eegfield)
                    EEGTMPTMP = std_getdataset(STUDY, ALLEEG, 'design', g.design, 'cell', cellInds(iDes));
                    if strcmpi(g.eegrmdata, 'on'), EEGTMPTMP.data = []; EEGTMPTMP.icaact = []; end
                    EEGTMPTMP.(g.eegfield) = tmpData;
                    EEGTMP(iDes) = EEGTMPTMP;
                elseif length(g.datafield) == 1
                    if ~ischar(tmpData.(g.datafield{1})), error('Field content cannot be a string'); end
                    data(iDes,:,:,:) = tmpData.(g.datafield{1});
                elseif isfield(tmpData, 'data') && isempty(g.datafield)
                    data(iDes,:,:,:) = tmpData.data;
                else
                    data(iDes) = tmpData;
                end
            end
            data = shiftdim(data,1);
            if ~isempty(g.eegfield)
                returndata{cInd,gInd} = EEGTMP;
            else
                returndata{cInd,gInd} = data;
            end
        end
    end
