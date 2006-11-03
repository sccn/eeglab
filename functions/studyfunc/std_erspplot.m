% std_erspplot() - plot STUDY cluster ERSPs. Displays either mean cluster ERSPs, 
%                  or else all cluster component ERSPs plus the mean cluster 
%                  ERSP in one figure per condition. The ERSPs can be plotted 
%                  only if component ERSPs were computed and saved in the 
%                  EEG datasets in the STUDY. These may either be computed 
%                  during pre-clustering using the gui-based function 
%                  pop_preclust(), or via the equivalent commandline functions 
%                  eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%           >> [STUDY] = std_erspplot(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
% Optional inputs:
%   'clusters' - [numeric vector|'all'??]  cluster numbers to plot.
%                Else?? 'all' -> plot all clusters in STUDY 
%                {default: 'all'}.
%   'comps'    - [numeric vector|'all'??]  cluster components to plot.
%                Else?? 'all' -> plot all cluster components 
%                {default: 'all'}.
%   'channels' - [numeric vector]  channels to plot. {default: all??}
%   'mode'     - ['together'|'apart'] plotting mode. In 'together' 
%                mode, the average ERSPs of the requested clusters|channels  
%                are plotted in the same figure - one per condition ??and 
%                group??. In 'apart' mode, component ERSPs for each
%                cluster (or channel) are plotted in a separate 
%                figure (per condition ""and group??) along with their mean 
%                ERSP. Note that for clusters, this option is irrelevant if 
%                component indices are provided as input (via 'comps' above) 
%                {default: 'together'??} 
%   'figure'  - ['on'|'off'] 'on' -> plot on a new figure; 
%                'off' -> plot on current figure. 
%                Note: 'off' is optional for one cluster in 'together' 
%                mode. This is useful for incomporating a cluster ERSP 
%                into a complex figure. In case of multiple conditions, 
%                only the first condition is displayed, but clicking on 
%                the figure will open a new figure with all conditions 
%                plotted separately {default: 'on'} 
% Output:
%   STUDY     - the input STUDY set structure with the plotted cluster 
%               mean ERSPs added to allow quick replotting 
% Example:
%        >> [STUDY] = std_erspplot(STUDY,ALLEEG, 'clusters', 'all', ...
%                                       'mode', 'together');
%           % Plot the mean ERSPs of all clusters in STUDY together 
%           % on the same figure. 
%
% See also: pop_clustedit(), pop_preclust(), eeg_createdata(), eeg_preclust(), pop_clustedit()
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
% Revision 1.29  2006/10/03 21:54:31  scott
% help msg edits -- ?? remain    -sm
%
% Revision 1.28  2006/10/03 18:31:49  scott
% help msg edits. ARNO - SEE ??  -sm
%
% Revision 1.27  2006/10/02 22:29:21  scott
% minor help edits
%
% Revision 1.26  2006/09/12 18:52:01  arno
% reprogram from scratch (statistics...), backward compatible
%
                            
function [STUDY allersp alltimes ] = std_erspstatplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erspstatplot;
    return;
end;

STUDY = pop_erspparams(STUDY, 'default');

opt = finputcheck( varargin, { 'channels'    'cell'    []              {};
                               'caxis'       'real'    []              [];
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'itc' 'ersp' } 'ersp';
                               'mode'        'string'  []              '';
                               'plottf'      'real'    []              [];
                               'timerange'   'real'    []              STUDY.etc.erspparams.timerange;
                               'freqrange'   'real'    []              STUDY.etc.erspparams.freqrange;
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       'integer' []              []; % for backward compatibility
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'subject'     'string'  []              '';
                               'statmode'    'string'  { 'subjects' 'common' 'trials' } STUDY.etc.erspparams.statmode}, 'std_erspstatplot');
if isstr(opt), error(opt); end;

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if ~isempty(opt.comps), 
    opt.subject = STUDY.datasetinfo( STUDY.cluster(opt.clusters).sets(1,opt.comps)).subject;
end;

if ~isempty(opt.subject), groupstats = 'off'; disp('No group statistics for single subject');
else                      groupstats = STUDY.etc.erspparams.groupstats;
end;
if ~isempty(opt.subject), condstats = 'off'; disp('No condition statistics for single subject');
else                      condstats = STUDY.etc.erspparams.condstats;
end;
plotcurveopt = { ...
   'ersplim',    eval( [ 'STUDY.etc.erspparams.' opt.datatype 'lim' ]), ...
   'threshold',  STUDY.etc.erspparams.threshold, ...
   'maskdata',   STUDY.etc.erspparams.maskdata, ...
   'groupstats',  groupstats, ...
   'condstats',   condstats, ...
   'statistics', STUDY.etc.erspparams.statistics };

if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', opt.datatype, 'timerange', opt.timerange, 'statmode', opt.statmode);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', opt.datatype, 'timerange', opt.timerange, 'statmode', opt.statmode);
end;

% plot single scalp map
% ---------------------
if ~isempty(opt.plotfreq)
    allersp = cell(size(STUDY.changrp(1).erspdata));
    for ind = 1:length(STUDY.changrp(1).erspdata(:))
        allersp{ind} = zeros([ size(STUDY.changrp(1).erspdata{1}) length(opt.channels)]);
        for index = 1:length(allinds)
            if ~isempty(opt.channels)
                allersp{ind}(:,:,index)  = STUDY.changrp(allinds(index)).erspdata{ind};
                allfreqs                 = STUDY.changrp(allinds(index)).erspfreqs;
                alltimes                 = STUDY.changrp(allinds(index)).ersptimes;
            else
                allersp{ind}(:,:,index)  = STUDY.cluster(allinds(index)).erspdata{ind};
                allfreqs                 = STUDY.cluster(allinds(index)).erspfreqs;
                alltimes                 = STUDY.changrp(allinds(index)).ersptimes;
            end;
        end;
        allersp{ind} = permute(allersp{ind}, [1 3 2]);
    end;
    %erspbase(:,2) = [];
    %erspbase(:,1) = [];
    
    % select individual subject
    % -------------------------
    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(allersp,1)
            for g = 1:size(allersp,2)
                allersp{c,g} = allersp{c,g}(:,:,subjind);
            end;
        end;
    end;

    [pgroup pcond pinter] = std_plot({ allfreqs alltimes }, allersp, 'condnames', STUDY.condition, 'subject', opt.subject, 'legend', opt.legend, ...
                                      'datatype', opt.datatype,'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottf, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, 'topovals', opt.plottf, plotcurveopt{:});
    return;
end;

opt.legend = 'off';
if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);

for index = 1:length(allinds)

    if length(allinds) > 1, subplot(nr,nc,index); end;
    if ~isempty(opt.channels)
        eval( [ 'allersp  = STUDY.changrp(allinds(index)).' opt.datatype 'data;' ]);
        eval( [ 'alltimes = STUDY.changrp(allinds(index)).' opt.datatype 'times;' ]);
        eval( [ 'allfreqs = STUDY.changrp(allinds(index)).' opt.datatype 'freqs;' ]);
    else
        eval( [ 'allersp  = STUDY.cluster(allinds(index)).' opt.datatype 'data;' ]);
        eval( [ 'alltimes = STUDY.cluster(allinds(index)).' opt.datatype 'times;' ]);
        eval( [ 'allfreqs = STUDY.cluster(allinds(index)).' opt.datatype 'freqs;' ]);
    end;

    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(allersp,1)
            for g = 1:size(allersp,2)
                allersp{c,g} = allersp{c,g}(:,subjind);
            end;
        end;
    end;

    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plot({ allfreqs alltimes }, allersp, 'condnames', STUDY.condition, 'subject', opt.subject, 'legend', opt.legend, ...
                                      'datatype', opt.datatype,'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottf, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
