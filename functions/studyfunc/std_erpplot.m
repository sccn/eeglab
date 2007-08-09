% std_erpplot() - Command line function to plot STUDY cluster component ERPs. Either 
%                 displays grand mean ERPs for all requested clusters in the same figure, 
%                 with ERPs for different conditions (if any) plotted in different colors. 
%                 Else, displays ERP for each specified cluster in separate figures 
%                 (per condition), each containing the cluster component ERPs plus 
%                 the grand mean cluster ERP (in bold). ERPs can be plotted only if 
%                 component ERPs were computed and saved in the STUDY EEG datasets. 
%                 These can be computed during pre-clustering using the gui-based 
%                 function pop_preclust() or the equivalent command line functions 
%                 eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
%                 and std_propplot().
% Usage:    
%   >> [STUDY] = std_erpplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY erpdata erptimes pgroup pcond pinter] = std_erpplot(STUDY, ALLEEG, ...);  
%
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets included 
%                in the STUDY. A STUDY set ALLEEG is typically created by load_ALLEEG().  
%
% Optional inputs for component plotting:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                ERPs of the requested clusters are plotted in the same figure, 
%                with ERPs for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, ERPS for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component ERPs plus the
%                average cluster ERP in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERP is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot ERP of all subjects.
%
% Other optional inputs:
%   'plotmode'  - ['normal'|'condensed'] 'normal'  -> plot in a new figure; 
%                 'condensed' -> plot all curves in the current figure in a 
%                 condensed fashion {default: 'normal'}
%   'key','val' - All optional inputs to pop_erpparams() are also accepted here
%                 to plot subset of time, statistics etc. The values used by default
%                 are the ones set using pop_erpparams() and stored in the
%                 STUDY structure.
%
% Outputs:
%   STUDY      - the input STUDY set structure with plotted cluster mean
%                ERPs data to allow quick replotting 
%   erpdata    - [cell] ERP data for each condition, group and subjects.
%                size of cell array is [nconds x ngroups]. Size of each element
%                is [times x subjects] for data channels or [times x components]
%                for component clusters. This array may be gicen as input 
%                directly to the statcond() function or std_stats function
%                to compute statistics.
%   erptimes   - [array] ERP time point latencies.
%   pgroup     - [array or cell] p-values group statistics. Output of the 
%                statcond() function.
%   pcond      - [array or cell] condition statistics. Output of the statcond() 
%                function.
%   pinter     - [array or cell] groups x conditions statistics. Output of
%                statcond() function.
%
%   Example:
%            >> [STUDY] = std_erpplot(STUDY,ALLEEG, 'clusters', 2, 'comps', 'all');
%               % Plot cluster-2 component ERPs plus the mean ERP in bold. 
%
%  See also  pop_clustedit(), pop_preclust(), eeg_createdata(), eeg_preclust(). std_propplot()
%
% Authors: Arnaud Delorme, CERCO, August, 2006-

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
% Revision 1.53  2007/08/09 18:29:18  arno
% add warning message and disable statistics for multiple cluster plotting
%
% Revision 1.52  2007/08/07 23:00:32  arno
% fix message
%
% Revision 1.51  2007/08/07 22:58:21  arno
% debug last change
%
% Revision 1.50  2007/08/07 22:55:27  arno
% do not reload the data to not reinvert polarities
%
% Revision 1.49  2007/06/25 04:33:44  toby
% altered multiple dipole plot windows to indicate number of components and subjects
%
% Revision 1.48  2007/05/25 03:23:51  toby
% oops, mispelled the correct index variable
%
% Revision 1.47  2007/05/25 03:18:26  toby
% polarity of component assigned by wrong index, causing crashes and incorrect polarity switching.
%
% Revision 1.46  2007/05/01 21:17:07  arno
% clarify help message
%
% Revision 1.45  2007/04/30 20:11:34  arno
% update header and code
%
% Revision 1.44  2007/04/28 00:28:21  arno
% backward compatibility
%
% Revision 1.43  2007/04/06 02:20:54  arno
% right unit
%
% Revision 1.42  2007/04/06 02:15:49  arno
% fix figure creation
%
% Revision 1.41  2007/04/05 22:01:21  arno
% condensed mode
%
% Revision 1.40  2007/03/20 03:30:59  arno
% fixed reading ERPs if icatopo absent
%
% Revision 1.39  2007/03/17 21:10:57  arno
% Matlab 6.5 compatibility
%
% Revision 1.38  2007/03/14 03:13:28  arno
% ERP polarity inversion
%
% Revision 1.37  2007/03/14 01:15:37  arno
% plot condensed mode
%
% Revision 1.36  2007/03/14 00:56:19  arno
% plotting ERP for multiple components
%
% Revision 1.35  2007/02/28 12:03:41  arno
% output statistics and documentation
%
% Revision 1.33  2007/01/26 18:04:33  arno
% reprogrammed from scractch (again)
%
% Revision 1.32  2006/11/23 00:26:16  arno
% add subject name
%
% Revision 1.31  2006/11/22 20:06:18  arno
% filter option
%
% Revision 1.30  2006/11/22 19:43:01  arno
% cannot select subject and component at the same time
%
% Revision 1.29  2006/11/09 22:44:06  arno
% copying changes to std_specplot
%
% Revision 1.28  2006/11/02 21:48:44  arno
% header
%
% Revision 1.27  2006/10/03 21:48:48  scott
% edit help msg -- some ?? remain   -sm
% ,
%
% Revision 1.26  2006/10/03 18:25:34  scott
% help msg edits ARNO - SEE ??  -sm
%
% Revision 1.25  2006/10/02 20:26:18  scott
% plotcond -> plotconditions
%
% Revision 1.24  2006/10/02 17:22:47  scott
% changed 'plotgroup' to 'plotgroups' ala change in std_erpparams.m
%
% Revision 1.23  2006/09/12 18:50:03  arno
% reprogram from scratch (statistics...), backward compatible
%
                            
function [STUDY, erpdata, alltimes, pgroup, pcond, pinter] = std_erpplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erpplot;
    return;
end;

STUDY = pop_erpparams(STUDY, 'default');

opt = finputcheck( varargin, { 'topotime'    'real'    [] STUDY.etc.erpparams.topotime;
                               'filter'      'real'    [] STUDY.etc.erpparams.filter;
                               'timerange'   'real'    [] STUDY.etc.erpparams.timerange;
                               'ylim'        'real'    [] STUDY.etc.erpparams.ylim;
                               'statistics'  'string'  [] STUDY.etc.erpparams.statistics;
                               'groupstats'  'string'  [] STUDY.etc.erpparams.groupstats;
                               'condstats'   'string'  [] STUDY.etc.erpparams.condstats;
                               'plotgroups'  'string'  [] STUDY.etc.erpparams.plotgroups;
                               'plotconditions' 'string'  [] STUDY.etc.erpparams.plotconditions;
                               'threshold'   'real'    [] STUDY.etc.erpparams.threshold;
                               'naccu'       'integer' [] STUDY.etc.erpparams.naccu;
                               'channels'    'cell'    []              {};
                               'caxis'       'real'    []              [];
                               'clusters'    'integer' []              [];
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       { 'string' 'integer' } [] []; % for backward compatibility
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'subject'     'string'  []              '';
                               'statmode'    'string'  { 'subjects' 'common' 'trials' } 'subjects'}, 'std_erpplot');
if isstr(opt), error(opt); end;
if isstr(opt.comps), opt.comps = []; opt.plotsubjects = 'on'; end;

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if ~isempty(opt.subject), opt.groupstats = 'off'; disp('No group statistics for single subject'); end;
if ~isempty(opt.subject), opt.condstats = 'off'; disp('No condition statistics for single subject'); end;
plotcurveopt = { ...
   'ylim',           opt.ylim, ...
   'filter',         opt.filter, ...
   'threshold',      opt.threshold, ...
   'plotgroups',     opt.plotgroups, ...
   'plotconditions', opt.plotconditions, ...
   'statistics',     opt.statistics };

if ~isnan(opt.topotime) & length(opt.channels) < 5
    warndlg2(strvcat('ERP parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;

% read data from disk
% -------------------
if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'erp', 'timerange', opt.timerange);
else 
     % invert polarity of ERPs
     filename = fullfile( ALLEEG(1).filepath, ALLEEG(1).filename(1:end-3));
     if exist([filename 'icatopo']) ~= 0
         for cls = opt.clusters
             reload = 0;
             if ~isfield(STUDY.cluster(opt.clusters(1)), 'erpdata'), reload = 1;
             elseif isempty(STUDY.cluster(opt.clusters(1)).erpdata), reload = 1;
             end;
             if reload
                 STUDY = std_readdata(STUDY, ALLEEG, 'clusters', cls, 'infotype', 'erp', 'timerange', opt.timerange);
                 STUDY = std_readtopoclust(STUDY, ALLEEG, cls);
                 if size(STUDY.cluster(opt.clusters(1)).erpdata,2) > 1 & cls == opt.clusters(1)
                     disp('WARNING: component polarity inversion for ERP not implemented if more than 1 group');
                     disp('WARNING: ERPs may not have the correct polarity');
                 elseif cls == opt.clusters(1)
                     disp('Inverting ERP component polarities based on scalp map polarities');
                 end    
                 clust = STUDY.cluster(cls);
                 for index = 1:length(clust.erpdata)
                     for comps = 1:size(clust.erpdata{index},2)
                         clust.erpdata{index}(:,comps) = clust.erpdata{index}(:,comps)*clust.topopol(comps);
                     end;
                 end;
                 STUDY.cluster(cls) = clust;
             end;
         end;
     end;
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'erp', 'timerange', opt.timerange);
end;

if strcmpi(opt.plotmode, 'condensed') | ...
   (length(allinds) > 1 & isempty(opt.channels))
    plotcurveopt = { plotcurveopt{:} 'figure' 'off' }; 
end;

% channel plotting
% ----------------
if ~isempty(opt.channels)
    structdat = STUDY.changrp;
    erpdata = cell(size(structdat(allinds(1)).erpdata));
    for ind =  1:length(structdat(allinds(1)).erpdata(:))
        erpdata{ind} = zeros([ size(structdat(allinds(1)).erpdata{1}) length(allinds)]);
        for index = 1:length(allinds)
            erpdata{ind}(:,:,index)  = structdat(allinds(index)).erpdata{ind};
            alltimes                 = structdat(allinds(index)).erptimes;
            compinds                 = structdat(allinds(index)).allinds;
            setinds                  = structdat(allinds(index)).setinds;
        end;
        erpdata{ind} = squeeze(permute(erpdata{ind}, [1 3 2])); % time elec subjects
    end;

    % select specific subject or component then plot
    % ----------------------------------------------
    if ~isempty(opt.subject), erpdata = std_selsubject(erpdata, opt.subject, setinds, { STUDY.datasetinfo(:).subject }, length(STUDY.subject)); end;
    
    % select specific time    
    % --------------------
    if ~isempty(opt.topotime)
        [tmp ti1] = min(abs(alltimes-opt.topotime(1)));
        [tmp ti2] = min(abs(alltimes-opt.topotime(end)));
        for index = 1:length(erpdata(:))
            erpdata{index} = mean(erpdata{index}(ti1:ti2,:,:),1);
        end;
        if opt.topotime(1) == opt.topotime(end), titlestr = [ num2str(opt.topotime(1)) ' ms'];
        else                                     titlestr = [ num2str(opt.topotime(1)) '-' num2str(opt.topotime(2)) ' ms'];
        end;
    end;
    
    % compute statistics and plot
    % ---------------------------
    [pcond pgroup pinter] = std_stat(erpdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, ...
                                         'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold);
    locs = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locs(std_chaninds(STUDY, opt.channels));
    if ~isempty(opt.topotime)
        std_chantopo(erpdata, 'condnames', STUDY.condition, 'plottopo', fastif(length(allinds)==1, 'off', 'on'), ...
                                      'datatype', 'spec', 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'unitx', '\muV', ...
                                      'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, 'topovals', titlestr, plotcurveopt{:});
    else
        std_plotcurve(alltimes, erpdata, 'condnames', STUDY.condition, 'plottopo', fastif(length(allinds)==1, 'off', 'on'), ...
                                      'datatype', 'spec', 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'unitx', 'ms', ...
                                      'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    end;
else 
    % plot component
    % --------------
    opt.legend = 'off';
    if length(allinds) > 1, figure('color', 'w'); end;
    nc = ceil(sqrt(length(allinds)));
    nr = ceil(length(allinds)/nc);
    comp_names = {};

    if length(opt.clusters) > 1 & isnan(opt.threshold) & ...
            ( strcmpi(opt.condstats, 'on') | strcmpi(opt.groupstats, 'on'))
        opt.condstats = 'off'; opt.groupstats = 'off'; 
        disp('Statistics disabled for plotting multiple clusters unless a threshold is set');
    end;

    for index = 1:length(allinds)

        if length(allinds) > 1, subplot(nr,nc,index); end;
        erpdata  = STUDY.cluster(allinds(index)).erpdata;
        alltimes = STUDY.cluster(allinds(index)).erptimes;
        compinds = STUDY.cluster(allinds(index)).allinds;
        setinds  = STUDY.cluster(allinds(index)).setinds;

        % plot specific component
        % -----------------------
        [erpdata opt.subject comp_names] = std_selcomp(STUDY, erpdata, allinds(index), setinds, compinds, opt.comps);
        [pcond pgroup pinter] = std_stat(erpdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, ...
                                         'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold);
            
        if index == length(allinds), opt.legend = 'on'; end;
            std_plotcurve(alltimes, erpdata, 'condnames', STUDY.condition, 'legend', opt.legend, 'subject', opt.subject, ...
                                          'compinds', comp_names, 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', ...
                          opt.topotime, 'unitx', 'ms',  'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
        if length(allinds) > 1, 
            if isempty(opt.channels), %title(sprintf('Cluster %d', allinds(index))); 
                title([ STUDY.cluster(allinds(index)).name ' (' num2str(length(STUDY.cluster(allinds(index)).comps)),...
                        ' ICs, '  num2str(length(unique(STUDY.cluster(allinds(index)).sets(1,:)))) ' Ss)' ]);
            else                      title(sprintf('%s', opt.channels{index}));  
            end;
        end;
    end;
end;
