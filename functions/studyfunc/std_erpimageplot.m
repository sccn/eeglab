% std_erpimageplot() - Commandline function to plot cluster ERPimage or channel erpimage.
%
% Usage:
%   >> [STUDY] = std_erpimageplot(STUDY, ALLEEG, key1, val1, key2, val2);
%   >> [STUDY data times freqs pgroup pcond pinter] = ...
%                std_erpimageplot(STUDY, ALLEEG ...);
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets in the STUDY.
%                Note: ALLEEG for a STUDY set is typically created using load_ALLEEG().
%
% Additional help:
% Inputs and output of this function are strictly identical to the std_erspplot().
% See the help message of this function for more information.
%
% See also: std_erspplot()
%
% Authors: Arnaud Delorme, UCSD/CERCO, August, 2011-

% Copyright (C) Arnaud Delorme, arno@ucsd.edu
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

function [STUDY, allitc, alltimes, allfreqs, pgroup, pcond, pinter] = std_erpimageplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erpimageplot;
    return;
end;

STUDY = pop_erpimparams(STUDY, 'default');
if strcmpi(STUDY.etc.erpimparams.concatenate, 'off')

    % use std_erspplot for stats when concatenate is off
    [STUDY allitc alltimes allfreqs pgroup pcond pinter ] = std_erspplot(STUDY, ALLEEG, 'datatype', 'erpim', varargin{:});

else
    params = STUDY.etc.erpimparams;
    [ opt moreparams ] = finputcheck( varargin, { ...
        'design'      'integer' [] STUDY.currentdesign;
        'topotime'    'real'    [] params.topotime;
        'topotrial'  'real'     [] params.topotrial;
        'timerange'   'real'    [] params.timerange;
        'trialrange'  'real'    [] params.trialrange;
        'colorlimits' 'real'    [] params.colorlimits; % ERPimage
        'statistics'  'string'  [] params.statistics;
        'groupstats'  'string'  [] params.groupstats;
        'condstats'   'string'  [] params.condstats;
        'threshold'   'real'    [] params.threshold;
        'naccu'       'integer' [] params.naccu;
        'mcorrect'    'string'  [] params.mcorrect;
        'erpimageopt' 'cell'    [] params.erpimageopt;
        'concatenate'  'string' { 'on','off' }  params.concatenate;
        'channels'    'cell'    []              {};
        'clusters'    'integer' []              [];
        'comps'       {'integer','string'}  []              []; % for backward compatibility
        'plotsubjects' 'string' { 'on','off' }  'off';
        'plotmode'    'string' { 'normal','condensed','none' }  'normal';
        'subject'     'string'  []              '' }, 'std_erpimageplot', 'ignore');
    
    if ~isempty(opt.topotime) && ~isempty(opt.topotrial)
        error('Cannot plot topography when ERP-image is in trial concatenation mode');
    end;
    if ~isempty(opt.trialrange)
        error('Cannot select trial range when ERP-image is in trial concatenation mode');
    end;
    if strcmpi(opt.groupstats, 'on') || strcmpi(opt.condstats, 'on') 
        disp('Warning: cannot perform statistics when ERP-image is in trial concatenation mode');
    end;
    
    % options
    if ~isempty(opt.colorlimits), options = { 'caxis' opt.colorlimits opt.erpimageopt{:} };
    else                          options = { 'cbar' 'on' opt.erpimageopt{:}  };
    end;
    
    if ~isempty(opt.channels)
        [STUDY allerpimage alltimes alltrials tmp events] = std_readersp(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'erpim', 'subject', opt.subject, ...
            'concatenate', 'on', 'timerange', opt.timerange, 'design', opt.design);
        
        % get figure title
        % ----------------
        locs = eeg_mergelocs(ALLEEG.chanlocs);
        locs = locs(std_chaninds(STUDY, opt.channels));
        allconditions = STUDY.design(opt.design).variable(1).value;
        allgroups     = STUDY.design(opt.design).variable(2).value;
        alltitles = std_figtitle('condnames', allconditions, 'cond2names', allgroups, 'chanlabels', { locs.labels }, ...
            'subject', opt.subject, 'valsunit', 'ms', 'datatype', 'ERPIM');
        
        figure;
        for iCond = 1:length(allconditions)
            for iGroup = 1:length(allgroups)
                tmpevents = events{iCond, iGroup};
                if isempty(tmpevents), tmpevents = zeros(1, size(allerpimage{iCond, iGroup},2)); end;
                subplot(length(allconditions), length(allgroups), (iCond-1)*length(allgroups) + iGroup);
                
                % use color scale for last plot
                if ~isempty(opt.colorlimits) && iCond == length(allconditions) && iGroup == length(allgroups)
                    options = { options{:} 'cbar' 'on' };
                end;
                erpimage(allerpimage{iCond, iGroup}, tmpevents, alltimes, alltitles{iCond, iGroup}, params.smoothing, params.nlines, options{:});
            end;
        end;
        
    else
        for cInd = 1:length(opt.clusters)
            [STUDY allerpimage alltimes alltrials tmp events] = std_readersp(STUDY, ALLEEG, 'clusters', opt.clusters(cInd), 'infotype', 'erpim', 'subject', opt.subject, ...
                'concatenate', 'on', 'timerange', opt.timerange, 'design', opt.design);
            % get figure title
            % ----------------
            locs = eeg_mergelocs(ALLEEG.chanlocs);
            locs = locs(std_chaninds(STUDY, opt.channels));
            allconditions = STUDY.design(opt.design).variable(1).value;
            allgroups     = STUDY.design(opt.design).variable(2).value;
            alltitles = std_figtitle('condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(cInd)).name, ...
                'subject', opt.subject, 'valsunit', 'ms', 'datatype', 'ERPIM');
            
            figure;
            for iCond = 1:length(allconditions)
                for iGroup = 1:length(allgroups)
                    tmpevents = events{iCond, iGroup};
                    if isempty(tmpevents), tmpevents = zeros(1, size(allerpimage{iCond, iGroup},2)); end;
                    subplot(length(allconditions), length(allgroups), (iCond-1)*length(allgroups) + iGroup);
                    
                    % use color scale for last plot
                    if ~isempty(opt.colorlimits) && iCond == length(allconditions) && iGroup == length(allgroups)
                        options = { options{:} 'cbar' 'on' };
                    end;
                    erpimage(allerpimage{iCond, iGroup}, tmpevents, alltimes, alltitles{iCond, iGroup}, params.smoothing, params.nlines, options{:});
                end;
            end;
        end;
    end;
end;
