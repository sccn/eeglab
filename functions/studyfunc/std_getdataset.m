% std_getdataset() - Constructs and returns EEG dataset from STUDY design.
%
% Usage:    
%     >> EEG = std_getdataset(STUDY, ALLEEG, 'key', 'val', ...);
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
%   'design'     - [numeric vector] STUDY design index. Default is to use
%                  the current design.
%   'rmcomps'    - [integer array] remove artifactual components (this entry
%                  is ignored when plotting components). This entry contains 
%                  the indices of the components to be removed. Default is none.
%   'rmclust'    - [integer array] which cluster(s) to remove from the data
%                  Default is none.
%   'interp'     - [struct] channel location structure containing electrode
%                  to interpolate ((this entry is ignored when plotting 
%                  components). Default is no interpolation.
%   'cell'       - [integer] index of the STUDY design cell to convert to
%                  an EEG dataset. Default is the first cell.
%   'onecomppercluster' - ['on'|'off'] when 'on' enforces one component per 
%                  cluster. Default is 'off'.
%   'interpcomponent' - ['on'|'off'] when 'on' interpolating component 
%                  scalp maps. Default is 'off'.
%   'cluster'    - [integer] which cluster(s). When this option is being
%                  used, only the component contained in the selected
%                  clusters are being loaded in the dataset.
%   'checkonly'  - ['on'|'off'] use in conjunction with the option above.
%                  When 'on', no dataset is returned. Default is 'off'.
%
% Outputs:
%   EEG          - EEG dataset structure corresponding to the selected
%                  STUDY design cell (element)
%   complist     - list of components selected
%
% Example to build the dataset corresponding to the first cell of the first
% design:
%   EEG = std_getdataset(STUDY, ALLEEG, 'design', 1, 'cell', 1);
%
% Example to check that all datasets in the design have exactly one
% component per cluster for cluster 2, 3, 4 and 5.
%   std_getdataset(STUDY, ALLEEG, 'design', 1, 'cell',
%   [1:length(STUDY.design(1).cell)], 'cluster', [2 3 4 5], 'checkonly', 'on');
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 12, 2004, arno@sccn.ucsd.edu
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

function [EEGOUT listcomp rmlistcomp] = std_getdataset(STUDY, ALLEEG, varargin);

if nargin < 2
    help std_getdataset;
    return;
end;

opt = finputcheck( varargin, { 'design'  'integer'   []    STUDY.currentdesign;
                               'interp'  'struct'    { }   struct([]);
                               'rmcomps' 'integer'	 []    [];
                               'cluster' 'integer'   []    [];
                               'rmclust' 'integer'   []    [];
                               'onecomppercluster' 'string'   {'on' 'off'}    'off';
                               'interpcomponent'   'string'   {'on' 'off'}    'off';
                               'checkonly' 'string'   {'on' 'off'}    'off';
                               'cell'    'integer'   []    1                    }, 'std_getdataset');
%                               'mode'   'string'   { 'channels' 'components' }  'channels';
if isstr(opt), error(opt); end;

if length(opt.cell) > 1
    % recursive call if more than one dataset
    % ---------------------------------------
    indcell = strmatch('cell', varargin(1:2:end));
    for index = 1:length(opt.cell)
        varargin{2*indcell} = index;
        EEGOUT(index) = std_getdataset(STUDY, ALLEEG, varargin{:});
    end;
else
    mycell = STUDY.design(opt.design).cell(opt.cell);
    
    % find components in non-artifactual cluster
    if ~isempty(opt.cluster)
        listcomp = [];
        for index = 1:length(opt.cluster)
            clsset  = STUDY.cluster(opt.cluster(index)).sets;
            clscomp = STUDY.cluster(opt.cluster(index)).comps;
            
            [indrow indcomp] = find(clsset == mycell.dataset(1));
            if length(indcomp) ~= 1 && strcmpi(opt.onecomppercluster, 'on')
                error(sprintf('Dataset %d must have exactly 1 components in cluster %d', mycell.dataset(1), opt.cluster(index)));
            end;
            listcomp = [listcomp clscomp(indcomp')];
        end;
    end;
    
    % find components in artifactual clusters
    if ~isempty(opt.rmclust)
        rmlistcomp = [];
        for index = 1:length(opt.rmclust)
            clsset  = STUDY.cluster(opt.rmclust(index)).sets;
            clscomp = STUDY.cluster(opt.rmclust(index)).comps;
            [indrow indcomp] = find(clsset == mycell.dataset(1));
            rmlistcomp = [rmlistcomp clscomp(indcomp')];
        end;
        if ~isempty(opt.rmcomps)
            disp('Both ''rmclust'' and ''rmcomps'' are being set. Artifact components will be merged');
        end;
        opt.rmcomps = [ opt.rmcomps rmlistcomp ];
    end;  
    rmlistcomp  = opt.rmcomps;
    if strcmpi(opt.checkonly, 'on'), EEGOUT = 0; return; end;
    
    % get data
    EEG = ALLEEG(mycell.dataset);
    EEGOUT = EEG(1);
    EEGOUT.data = eeg_getdatact(EEG, 'channel', [], 'trialindices', mycell.trials, 'rmcomps', opt.rmcomps, 'interp', opt.interp);
    EEGOUT.trials   = size(EEGOUT.data,3);
    EEGOUT.nbchan   = size(EEGOUT.data,1);
    if ~isempty(opt.interp)
        EEGOUT.chanlocs = opt.interp;
    end;
    EEGOUT.event  = [];
    EEGOUT.epoch  = [];
    EEGOUT.filename = mycell.filebase;
    EEGOUT.condition = mycell.value{1};
    EEGOUT.group     = mycell.value{2};
    EEGOUT.subject   = mycell.case;
    if ~isempty(opt.cluster)
        if ~isempty(opt.interp) && strcmpi(opt.interpcomponent, 'on')
            TMPEEG = EEGOUT;
            TMPEEG.chanlocs = EEG(1).chanlocs;
            TMPEEG.data = EEG(1).icawinv;
            TMPEEG.nbchan = size(TMPEEG.data,1);
            TMPEEG.pnts   = size(TMPEEG.data,2);
            TMPEEG.trials = 1;
            TMPEEG = eeg_interp(TMPEEG, opt.interp, 'spherical');
            EEGOUT.icawinv = TMPEEG.data(:, listcomp);
        else
            EEGOUT.icawinv = EEGOUT.icawinv(:, listcomp);
        end;
        EEGOUT.icaact         = eeg_getdatact(EEG, 'component', listcomp, 'trialindices', mycell.trials );
        EEGOUT.icaweights     = EEGOUT.icaweights(listcomp,:);
        EEGOUT.etc.clustid    = { STUDY.cluster(opt.cluster).name };  % name of each cluster
        EEGOUT.etc.clustcmpid = listcomp;                             % index of each component in original ICA matrix
    else
        EEGOUT = eeg_checkset(EEGOUT);
    end;
end;
