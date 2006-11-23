% std_specplot() - plot STUDY component cluster spectra, either mean spectra 
%                  for all requested clusters in the same figure, with spectra 
%                  for different conditions (if any) plotted in different colors, 
%                  or spectra for each specified cluster in a separate figure 
%                  for each condition,  showing the cluster component spectra plus 
%                  the mean cluster spectrum (in bold). The spectra can be 
%                  plotted only if component spectra have been computed and 
%                  saved with the EEG datasets in Matlab files "[datasetname].icaspec" 
%                  using pop_preclust() or std_preclust(). Called by pop_clustedit(). 
%                  Calls std_readspec() and internal function std_plotcompspec()
% Usage:    
%              >> [STUDY] = std_specplot(STUDY, ALLEEG, key1, val1, key2, val2, ...);  
% Inputs:
%   STUDY      - STUDY structure comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - vector of EEG dataset structures for the dataset(s) in the STUDY, 
%                typically created using load_ALLEEG().  
% Optional inputs:
%   'clusters' - [int vector|'all'] -> cluster numbers to plot.
%                'all' -> plot all clusters in STUDY {default: 'all'}.
%   'comps'    - [int vector|'all'] -> indices of cluster components to plot.
%                'all' -> plot all the components in the cluster {default: 'all'}.
%   'mode'     - ['together'|'apart'] plotting mode. In 'together' mode, the average 
%                spectra of the requested clusters are plotted in the same figure, 
%                with spectra for  different conditions ??and groups?? (if any) 
%                plotted in different colors. In 'apart' mode, spectra for each 
%                specified cluster are plotted in separate figures (per condition 
%                ??and group??), 
%                each containing the overplotted individual cluster component spectra 
%                plus the mean cluster spectrum in bold.
%                Note that this option is irrelevant when component indices are provided 
%                as input (via 'apart' above) {default: 'together'}. 
%   'figure'   - ['on'|'off'] in 'together' mode, 'on' plots in a new figure, 
%                while 'off' plots in the current figure {default: 'on'}
% Outputs:
%   STUDY      - the input STUDY set structure with the plotted cluster mean spectra
%                added?? to allow quick replotting.
%   Example:
%            >> [STUDY] = std_specplot(STUDY,ALLEEG, 'clusters', 2, 'mode', 'apart');
%               % Plot component spectra for STUDY cluster 2, plus the mean cluster 
%               % spectrum (in bold). 
%
%  See also  pop_clustedit(), pop_preclust() std_preclust(), pop_clustedit(), std_readspec()
%
% Authors: Arnaud Delorme, CERCO, August, 2006

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

% $Log: not supported by cvs2svn $
% Revision 1.31  2006/11/22 19:41:17  arno
% cannot select subject and component at the same time
%
% Revision 1.30  2006/11/09 22:38:03  arno
% fix component and subject selection
%
% Revision 1.29  2006/10/04 23:46:58  toby
% Bug fix courtesy Bas de Kruif
%
% Revision 1.28  2006/10/03 22:22:35  scott
% help msg
%
% Revision 1.27  2006/10/03 21:57:50  scott
% edit help msg  -- some ?? remain   -sm
%
% Revision 1.26  2006/10/03 18:39:25  scott
% help msg   ARNO - SEE ??    -sm
%
% Revision 1.25  2006/10/03 18:37:49  scott
% help msg edits.  ARNO - SEE ??   -sm
%
% Revision 1.24  2006/10/02 20:25:30  scott
% plotcond -> plotconditions
%
% Revision 1.23  2006/10/02 17:24:38  scott
% edited help msg for clarity, changed 'plotgroup' to 'plotgroups'
% ala change in std_specparams
%
% Revision 1.22  2006/10/02 11:43:00  arno
% allow plotting scalp maps
%
                            
function [STUDY, erspbase] = std_specplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_specplot;
    return;
end;

STUDY = pop_specparams(STUDY, 'default');

opt = finputcheck( varargin, { 'channels'    'cell'   []              {};
                               'caxis'       'real'   []              [];
                               'clusters'    'integer' []              [];
                               'plotfreq'    'real'   []              [];
                               'freqrange'   'real'   []              STUDY.etc.specparams.freqrange;
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       'integer' []              []; % for backward compatibility
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'subject'     'string' []              '';
                               'statmode'    'string' { 'individual' 'common' 'trials' } 'individual'}, 'std_specplot');

if isstr(opt), error(opt); end;

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'apart'), opt.plotsubjects = 'on'; end;
if ~isempty(opt.comps), 
    opt.subject = STUDY.datasetinfo( STUDY.cluster(opt.clusters).sets(1,opt.comps)).subject;
end;

if ~isempty(opt.subject), groupstats = 'off'; disp('No group statistics for single subject');
else                      groupstats = STUDY.etc.specparams.groupstats;
end;
if ~isempty(opt.subject), condstats = 'off'; disp('No condition statistics for single subject');
else                      condstats = STUDY.etc.specparams.condstats;
end;
plotcurveopt = { ...
   'ylim',       STUDY.etc.specparams.ylim, ...
   'threshold',  STUDY.etc.specparams.threshold, ...
   'groupstats',  groupstats, ...
   'condstats',   condstats, ...
   'plotgroups',  STUDY.etc.specparams.plotgroups, ...
   'plotconditions', STUDY.etc.specparams.plotconditions, ...
   'statistics', STUDY.etc.specparams.statistics };

if ~isempty(opt.plotfreq) & ~isempty(opt.channels)
    alllocs      = eeg_mergelocs(ALLEEG(:).chanlocs);
    opt.channels = { alllocs.labels };
end;
if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'spec', 'freqrange', opt.freqrange);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'spec', 'freqrange', opt.freqrange);
end;

% plot single scalp map
% ---------------------
if ~isempty(opt.plotfreq)
    erspbase = cell(size(STUDY.changrp(1).specdata));
    for ind = 1:length(STUDY.changrp(1).specdata(:))
        erspbase{ind} = zeros([ size(STUDY.changrp(1).specdata{1}) length(opt.channels)]);
        for index = 1:length(allinds)
            if ~isempty(opt.channels)
                erspbase{ind}(:,:,index) = STUDY.changrp(allinds(index)).specdata{ind};
                allfreqs                 = STUDY.changrp(allinds(index)).specfreqs;
            else
                erspbase{ind}(:,:,index) = STUDY.cluster(allinds(index)).specdata{ind};
                allfreqs                 = STUDY.cluster(allinds(index)).specfreqs;
            end;
        end;
        erspbase{ind} = permute(erspbase{ind}, [1 3 2]);
    end;
    %erspbase(:,2) = [];
    %erspbase(:,1) = [];
    
    % select individual subject
    % -------------------------
    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(erspbase,1)
            for g = 1:size(erspbase,2)
                erspbase{c,g} = erspbase{c,g}(:,:,subjind);
            end;
        end;
    end;

    [pgroup pcond pinter] = std_plot(allfreqs, erspbase, 'condnames', STUDY.condition, ...
                                      'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plotfreq, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    return;
end;

opt.legend = 'off';
if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);
comp_names = {};

for index = 1:length(allinds)

    if length(allinds) > 1, subplot(nr,nc,index); end;
    if ~isempty(opt.channels)
        erspbase = STUDY.changrp(allinds(index)).specdata;
        allfreqs = STUDY.changrp(allinds(index)).specfreqs;
        setinds  = STUDY.changrp(allinds(index)).setinds;
    else
        erspbase = STUDY.cluster(allinds(index)).specdata;
        allfreqs = STUDY.cluster(allinds(index)).specfreqs;
        compinds = STUDY.cluster(allinds(index)).allinds;
        setinds  = STUDY.cluster(allinds(index)).setinds;
    end;
    
    % plot specific subject
    % ---------------------
    if ~isempty(opt.subject) & isempty(opt.comps)
        for c = 1:size(erspbase,1)
            for g = 1:size(erspbase,2)
                for l=length(setinds{c,g}):-1:1
                    if ~strcmpi(opt.subject, STUDY.datasetinfo(setinds{c,g}(l)).subject)
                        erspbase{c,g}(:,l) = [];
                    end;
                end;
            end;
        end;
    end;
    
    % plot specific component
    % -----------------------
    if ~isempty(opt.comps)
        
        % find and select group
        % ---------------------
        sets   = STUDY.cluster(allinds(index)).sets(:,opt.comps);
        comps  = STUDY.cluster(allinds(index)).comps(opt.comps);
        grp    = STUDY.datasetinfo(sets(1)).group;
        grpind = strmatch( grp, STUDY.group );
        if isempty(grpind), grpind = 1; end;
        erspbase = erspbase(:,grpind);
            
        % find component
        % --------------
        for c = 1:size(erspbase,1)
            for ind = 1:length(compinds{1,grpind})
                if compinds{1,grpind}(ind) == comps & any(setinds{1,grpind}(ind) == sets)
                    erspbase{c} = erspbase{c}(:,ind);
                    comp_names{c,1} = comps;
                end;
            end;
        end;
        opt.subject = STUDY.datasetinfo(sets(1)).subject;
    end;

    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plot(allfreqs, erspbase, 'condnames', STUDY.condition, 'legend', opt.legend, 'subject', opt.subject, ...
                                      'compinds', comp_names, 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'plottopo', opt.plotfreq, 'unitx', 'Hz', ...
                                       'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
