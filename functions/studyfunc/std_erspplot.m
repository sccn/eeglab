% std_erspplot() - plot STUDY cluster ERSPs. Displays either mean cluster ERSPs, 
%                  or else all cluster component ERSPs plus the mean cluster 
%                  ERSP in one figure per condition. The ERSPs can be plotted 
%                  only if component ERSPs were computed and saved in the 
%                  EEG datasets in the STUDY. These may either be computed 
%                  during pre-clustering using the gui-based function 
%                  pop_preclust(), or via the equivalent commandline functions 
%                  eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%   >> [STUDY] = std_erspplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY erspdata ersptimes erspfreqs pgroup pcond pinter] = ...
%                std_erspplot(STUDY, ALLEEG ...);
%
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
% Optional inputs:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                ERSPs of the requested clusters are plotted in the same figure, 
%                with ERSPs for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, ERSP for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component ERSP plus the
%                average cluster ERSP in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERSP is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot ERSP of all subjects.
%
% Other optional inputs:
%   'figure'   - ['on'|'off'] 'on'  -> plot in a new figure; 
%                'off' -> plot in the current figure {default: 'on'}
%   'key','val' - All optional inputs to pop_specparams() are also accepted here
%                 to plot subset of time, statistics etc. The values used by default
%                 are the ones set using pop_specparams() and stored in the
%                 STUDY structure.
% Output:
%   STUDY      - the input STUDY set structure with the plotted cluster 
%                mean ERSPs added to allow quick replotting 
%   erspdata   - [cell] ERSP data for each condition, group and subjects.
%                size of cell array is [nconds x ngroups]. Size of each element
%                is [freqs x times x subjects] for data channels or 
%                [freqs x times x components] for component clusters. This 
%                array may be gicen as input  directly to the statcond() f
%                unction or std_stats() function to compute statistics.
%   ersptimes  - [array] ERSP time point latencies.
%   erspfreqs  - [array] ERSP point frequency values.
%   pgroup     - [array or cell] p-values group statistics. Output of the 
%                statcond() function.
%   pcond      - [array or cell] condition statistics. Output of the statcond() 
%                function.
%   pinter     - [array or cell] groups x conditions statistics. Output of
%                statcond() function.
%
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
% Revision 1.43  2007/02/28 12:04:12  arno
% output statistics and documentation
%
% Revision 1.42  2007/01/26 18:05:24  arno
% new handling of scalp plotting
%
% Revision 1.41  2006/11/23 00:28:32  arno
% add subject info to subject plot
%
% Revision 1.40  2006/11/22 23:01:07  arno
% array size for one subject
%
% Revision 1.39  2006/11/22 19:43:33  arno
% cannot select subject and component at the same time
%
% Revision 1.38  2006/11/15 21:29:39  arno
% subject index
%
% Revision 1.37  2006/11/15 01:59:23  arno
% dataset typo
% /
%
% Revision 1.36  2006/11/14 04:16:43  arno
% setinds for channels
%
% Revision 1.35  2006/11/09 22:04:35  arno
% fix cluster plotting
%
% Revision 1.34  2006/11/08 23:21:41  arno
% fixed plotting single subject
%
% Revision 1.33  2006/11/03 03:01:16  arno
% allow plotting specific time-freq point
%
% Revision 1.32  2006/11/03 02:11:22  arno
% same
%
% Revision 1.31  2006/11/03 02:09:47  arno
% same
%
% Revision 1.30  2006/11/03 02:08:47  arno
% allowing ploting single time-freq point
%
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
                            
function [STUDY, allersp, alltimes, allfreqs, pgroup, pcond, pinter] = std_erspplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erspstatplot;
    return;
end;

STUDY = pop_erspparams(STUDY, 'default');

opt = finputcheck( varargin, { 'topofreq'    'real'    [] STUDY.etc.specparams.topofreq;
                               'topotime'    'real'    [] STUDY.etc.specparams.topotime;
                               'freqrange'   'real'    [] STUDY.etc.specparams.freqrange;
                               'timerange'   'real'    [] STUDY.etc.specparams.timerange;
                               'ersplim'     'real'    [] STUDY.etc.specparams.ersplim;
                               'itclim'      'real'    [] STUDY.etc.specparams.itclim;
                               'subbaseline' 'string'  [] STUDY.etc.specparams.subbaseline;
                               'maskdata'    'string'  [] STUDY.etc.specparams.maskdata;
                               'statistics'  'string'  [] STUDY.etc.specparams.statistics;
                               'groupstats'  'string'  [] STUDY.etc.specparams.groupstats;
                               'condstats'   'string'  [] STUDY.etc.specparams.condstats;
                               'threshold'   'real'    [] STUDY.etc.specparams.threshold;
                               'naccu'       'integer' [] STUDY.etc.specparams.naccu;
                               'channels'    'cell'    []              {};
                               'caxis'       'real'    []              [];
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'itc' 'ersp' } 'ersp';
                               'mode'        'string'  []              '';
                               'plottf'      'real'    []              [];
                               'comps'       'integer' []              []; % for backward compatibility
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'subject'     'string'  []              '';
                               'statmode'    'string'  { 'subjects' 'common' 'trials' } STUDY.etc.erspparams.statmode}, 'std_erspstatplot');
if isstr(opt), error(opt); end;

% for backward compatibility
% --------------------------
if ~isnan(STUDY.etc.erspparams.topotime),
    opt.plottf = [ STUDY.etc.erspparams.topofreq STUDY.etc.erspparams.topotime ];
end;
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;

if ~isempty(opt.subject), opt.groupstats = 'off'; disp('No group statistics for single subject'); end;
if ~isempty(opt.subject), opt.condstats = 'off'; disp('No condition statistics for single subject'); end;

plotcurveopt = { ...
   'ersplim',     fastif(strcmpi(opt.datatype, 'ITC', opt.itclim, opt.ersplim), ...
   'threshold',   opt.threshold, ...
   'maskdata',    opt.maskdata, ...
   'groupstats',  opt.groupstats, ...
   'condstats',   opt.condstats, ...
   'statistics',  opt.statistics };
if ~isempty(opt.plottf) & length(opt.channels) < 5
    warndlg2(strvcat('ERSP/ITC parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;    

if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', opt.datatype, 'statmode', opt.statmode);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', opt.datatype, 'statmode', opt.statmode);
end;
opt.legend = 'off';

% plot single scalp map
% ---------------------
if ~isempty(opt.plottf)
    allersp = cell(size(STUDY.changrp(1).erspdata));
    for ind = 1:length(STUDY.changrp(1).erspdata(:))
        if size(STUDY.changrp(1).erspdata{1},3) == 1
             allersp{ind} = zeros([ size(STUDY.changrp(1).erspdata{1}) 1 length(opt.channels)]);
        else allersp{ind} = zeros([ size(STUDY.changrp(1).erspdata{1}) length(opt.channels)]);
        end;
        for index = 1:length(allinds)
            if ~isempty(opt.channels)
                eval( [ 'allersp{ind}(:,:,:,index)  = STUDY.changrp(allinds(index)).' opt.datatype 'data{ind};' ]);
                eval( [ 'alltimes                   = STUDY.changrp(allinds(index)).' opt.datatype 'times;' ]);
                eval( [ 'allfreqs                   = STUDY.changrp(allinds(index)).' opt.datatype 'freqs;' ]);
            else % not sure clusters are actually of any use here
                eval( [ 'allersp{ind}(:,:,:,index)  = STUDY.cluster(allinds(index)).' opt.datatype 'data{ind};' ]);
                eval( [ 'alltimes                   = STUDY.cluster(allinds(index)).' opt.datatype 'times;' ]);
                eval( [ 'allfreqs                   = STUDY.cluster(allinds(index)).' opt.datatype 'freqs;' ]);
            end;
        end;
        allersp{ind} = permute(allersp{ind}, [1 2 4 3]);
    end;
    %erspbase(:,2) = [];
    %erspbase(:,1) = [];
    
    % select individual subject
    % -------------------------
    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(allersp,1)
            for g = 1:size(allersp,2)
                allersp{c,g} = allersp{c,g}(:,:,:,subjind);
            end;
        end;
    end;
    
    locs = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locs(std_chaninds(STUDY, opt.channels));
    [pgroup pcond pinter] = std_plottf(alltimes, allfreqs, allersp, 'condnames', STUDY.condition, 'subject', opt.subject, 'legend', opt.legend, ...
                                      'datatype', opt.datatype,'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottf, 'unitx', 'Hz', ...
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, 'topovals', opt.plottf, plotcurveopt{:});
    return;
end;

if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);
comp_names = {};

for index = 1:length(allinds)

    if length(allinds) > 1, subplot(nr,nc,index); end;
    if ~isempty(opt.channels)
        eval( [ 'allersp  = STUDY.changrp(allinds(index)).' opt.datatype 'data;' ]);
        eval( [ 'alltimes = STUDY.changrp(allinds(index)).' opt.datatype 'times;' ]);
        eval( [ 'allfreqs = STUDY.changrp(allinds(index)).' opt.datatype 'freqs;' ]);
        setinds  = STUDY.changrp(allinds(index)).setinds;
    else
        eval( [ 'allersp  = STUDY.cluster(allinds(index)).' opt.datatype 'data;' ]);
        eval( [ 'alltimes = STUDY.cluster(allinds(index)).' opt.datatype 'times;' ]);
        eval( [ 'allfreqs = STUDY.cluster(allinds(index)).' opt.datatype 'freqs;' ]);
        compinds = STUDY.cluster(allinds(index)).allinds;
        setinds  = STUDY.cluster(allinds(index)).setinds;
    end;

    % plot specific subject
    % ---------------------
    if ~isempty(opt.subject) & isempty(opt.comps)
        for c = 1:size(allersp,1)
            for g = 1:size(allersp,2)
                for l=length(setinds{c,g}):-1:1
                    if ~strcmpi(opt.subject, STUDY.datasetinfo(setinds{c,g}(l)).subject)
                        allersp{c,g}(:,:,l) = [];
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
        allersp = allersp(:,grpind);
            
        % find component
        % --------------
        for c = 1:size(allersp,1)
            for ind = 1:length(compinds{1,grpind})
                if compinds{1,grpind}(ind) == comps & any(setinds{1,grpind}(ind) == sets)
                    allersp{c} = allersp{c}(:,:,ind);
                    comp_names{c,1} = comps;
                end;
            end;
        end;
        opt.subject = STUDY.datasetinfo(sets(1)).subject;
    end;
    
    % plot specific component
    % -----------------------
    
    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plottf(alltimes, allfreqs, allersp, 'condnames', STUDY.condition, 'subject', opt.subject, 'legend', opt.legend, ...
                                      'compinds', comp_names, 'datatype', opt.datatype,'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottf, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
