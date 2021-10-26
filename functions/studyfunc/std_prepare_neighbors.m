% std_prepare_neighbors() - prepare Fieldtrip channel neighbor structure.
%                           Only prepare the structure if necessary based
%                           on statistical options in STUDY.etc.statistics.
%                           Use the 'force' option to force preparing the
%                           matrix.
%
% Usage:
%   >> [STUDY neighbors] = std_prepare_neighbors( STUDY, ALLEEG, 'key', val)
%
% Inputs:
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%
% Optional inputs:
%   'force'        - ['on'|'off'] force generating the structure irrespective
%                    of the statistics options selected. Default is 'off'.
%   'channels'     - [cell] list of channels to include in the matrix
%   'method'       - [string] method for preparing. See ft_prepare_neighbors
%   'neighbordist' - [float] max distance. See ft_prepare_neighbors
%
% Note: other ft_prepare_neighbors fields such as template, layout may
%       also be used as optional keywords.
%
% Outputs:
%   STUDY       - an EEGLAB STUDY set of loaded EEG structures
%   neighbors   - Fieldtrip channel neighbour structure
%   limostruct  - structure compatible with LIMO (need to be saved in a file)
%
% Author: Arnaud Delorme, SCCN, UCSD, 2012-
%
% See also: statcondfieldtrip()

% Copyright (C) Arnaud Delorme
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

function [STUDY neighbors limostruct] = std_prepare_neighbors(STUDY, ALLEEG, varargin);

neighbors = [];
if nargin < 2
    return;
end

[opt addopts] = finputcheck( varargin, {  'force'    'string'  { 'on','off' }   'off';
                                          'channels' 'cell'    {}               {} }, 'std_stat', 'ignore');

if strcmpi(opt.force, 'on') || (strcmpi(STUDY.etc.statistics.fieldtrip.mcorrect, 'cluster') && ...
        strcmpi(STUDY.etc.statistics.mode, 'fieldtrip') && (strcmpi(STUDY.etc.statistics.groupstats, 'on') || strcmpi(STUDY.etc.statistics.condstats, 'on')))
    
    EEG = eeg_emptyset;
    EEG.chanlocs = eeg_mergelocs(ALLEEG.chanlocs);
    if isempty(EEG.chanlocs)
        disp('std_prepare_neighbors: cannot prepare channel neighbour structure because of empty channel structures');
        return;
    end
    
    if ~isempty(STUDY.etc.statistics.fieldtrip.channelneighbor) && isempty(addopts) && ...
        length(STUDY.etc.statistics.fieldtrip.channelneighbor) == length(EEG.chanlocs)
        
        disp('Using stored channel neighbour structure');
        neighbors = STUDY.etc.statistics.fieldtrip.channelneighbor;
        
    else

        if ~isempty(opt.channels)
            indChans = eeg_chaninds(EEG, opt.channels);
            EEG.chanlocs = EEG.chanlocs(indChans);
        end
        
        EEG.nbchan   = length(EEG.chanlocs);
        EEG.data     = zeros(EEG.nbchan,100,1);
        EEG.trials   = 1;
        EEG.pnts     = 100;
        EEG.xmin     = 0;
        EEG.srate    = 1;
        EEG.xmax     = 99;
        EEG = eeg_checkset(EEG);
        tmpcfg = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
        
        % call the function that find channel neighbors
        % ---------------------------------------------
        addparams = eval( [ '{' STUDY.etc.statistics.fieldtrip.channelneighborparam '}' ]);
        for index = 1:2:length(addparams)
            tmpcfg = setfield(tmpcfg, addparams{index}, addparams{index+1});
        end
        for index = 1:2:length(addopts)
            tmpcfg = setfield(tmpcfg, addopts{index}, addopts{index+1});
        end
        warning off;
        if isfield(EEG.chanlocs, 'theta') && ~isempty(EEG.chanlocs(1).theta)
            tmpcfg = rmfield(tmpcfg, 'cfg');
            tmpcfg2 = tmpcfg;
            tmpcfg  = rmfield(tmpcfg, 'label'); % first input must not be data
            tmpcfg2 = rmfield(tmpcfg2, 'method'); % second input must not be method
            if isfield(tmpcfg, 'trialinfo')
                tmpcfg  = rmfield(tmpcfg, 'trialinfo'); % first input must not be data
            end
            % tmpcfg = rmfield(tmpcfg, 'label');
            % --> removing label seems to make ft_prepare_neighbours to crash
            cfg.neighbors = ft_prepare_neighbours(tmpcfg, tmpcfg2);
            neighbors = cfg.neighbors;
        else
            neighbors = [];
        end
        warning on;
        
    end

    STUDY.etc.statistics.fieldtrip.channelneighbor = neighbors;
    
    if nargout > 2
        limostruct.expected_chanlocs = EEG.chanlocs;

        if ~isempty(neighbors)
            allLabels = { neighbors.label };
            limostruct.channeighbstructmat = zeros(length(EEG.chanlocs));
            for iN = 1:length(neighbors)
                if ~isequal(neighbors(iN).label, limostruct.expected_chanlocs(iN).labels)
                    error('Wrong label');
                else
                    [tmp posChan] = intersect( allLabels, neighbors(iN).neighblabel);
                    limostruct.channeighbstructmat(iN,posChan) = 1;
                    limostruct.channeighbstructmat(posChan,iN) = 1;
                end
            end
        else
            limostruct.channeighbstructmat = ones(length(EEG.chanlocs));
        end
    end
end
