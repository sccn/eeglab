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
%   'comps'    - [numeric vector|'all']  indices of the cluster components to plot.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERP is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%
% Other optional inputs:
%   'figure'   - ['on'|'off'] 'on'  -> plot in a new figure; 
%                'off' -> plot in the current figure {default: 'on'}
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
                            
function [STUDY erpdata alltimes pgroup pcond pinter] = std_erpplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erpplot;
    return;
end;

STUDY = pop_erpparams(STUDY, 'default');

opt = finputcheck( varargin, { 'channels'    'cell'    []              {};
                               'caxis'       'real'    []              [];
                               'clusters'    'integer' []              [];
                               'plottime'    'real'    []              [];
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       'integer' []              []; % for backward compatibility
                               'timerange'   'real'    []              STUDY.etc.erpparams.timerange;
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'subject'     'string'  []              '';
                               'statmode'    'string'  { 'subjects' 'common' 'trials' } 'subjects'}, 'std_erpplot');
if isstr(opt), error(opt); end;

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;

if isempty(opt.plottime), opt.plottime = STUDY.etc.erpparams.topotime; end;
if isnan(  opt.plottime), opt.plottime = []; end;
if ~isempty(opt.subject), groupstats = 'off'; disp('No group statistics for single subject');
else                      groupstats = STUDY.etc.erpparams.groupstats;
end;
if ~isempty(opt.subject), condstats = 'off'; disp('No condition statistics for single subject');
else                      condstats = STUDY.etc.erpparams.condstats;
end;
plotcurveopt = { ...
   'ylim',       STUDY.etc.erpparams.ylim, ...
   'filter',     STUDY.etc.erpparams.filter, ...
   'threshold',  STUDY.etc.erpparams.threshold, ...
   'plotgroups',  STUDY.etc.erpparams.plotgroups, ...
   'plotconditions',   STUDY.etc.erpparams.plotconditions, ...
   'statistics', STUDY.etc.erpparams.statistics };

if ~isempty(opt.plottime) & length(opt.channels) < 5
    warndlg2(strvcat('ERP parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;

% read data from disk
% -------------------
if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'erp', 'timerange', opt.timerange);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'erp', 'timerange', opt.timerange);
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
    if ~isempty(opt.plottime)
        [tmp ti1] = min(abs(alltimes-opt.plottime(1)));
        [tmp ti2] = min(abs(alltimes-opt.plottime(end)));
        for index = 1:length(erpdata(:))
            erpdata{index} = mean(erpdata{index}(ti1:ti2,:,:),1);
        end;
        if opt.plottime(1) == opt.plottime(end), titlestr = [ num2str(opt.plottime(1)) ' ms'];
        else                                     titlestr = [ num2str(opt.plottime(1)) '-' num2str(opt.plottime(2)) ' ms'];
        end;
    end;
    
    % compute statistics and plot
    % ---------------------------
    [pcond pgroup pinter] = std_stat(erpdata, STUDY.etc.erpparams, 'groupstats', groupstats, 'condstats', condstats);
    locs = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locs(std_chaninds(STUDY, opt.channels));
    if ~isempty(opt.plottime)
        std_chantopo(erpdata, 'condnames', STUDY.condition, 'plottopo', fastif(length(allinds)==1, 'off', 'on'), ...
                                      'datatype', 'spec', 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'unitx', '\muV', ...
                                      'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, 'topovals', titlestr, plotcurveopt{:});
    else
        std_plotcurve(alltimes, erpdata, 'condnames', STUDY.condition, 'plottopo', fastif(length(allinds)==1, 'off', 'on'), ...
                                      'datatype', 'spec', 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'unitx', '\muV', ...
                                      'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    end;
else 
    % plot component
    % --------------
    opt.legend = 'off';
    if length(allinds) > 1, figure('color', 'w'); opt.plotmode = 'condensed'; plotcurveopt = { plotcurveopt{:} 'figure' 'off' }; end;
    nc = ceil(sqrt(length(allinds)));
    nr = ceil(length(allinds)/nc);
    comp_names = {};

    for index = 1:length(allinds)

        if length(allinds) > 1, subplot(nr,nc,index); end;
        erpdata  = STUDY.cluster(allinds(index)).erpdata;
        alltimes = STUDY.cluster(allinds(index)).erptimes;
        compinds = STUDY.cluster(allinds(index)).allinds;
        setinds  = STUDY.cluster(allinds(index)).setinds;

        % plot specific component
        % -----------------------
        [erpdata opt.subject comp_names] = std_selcomp(STUDY, erpdata, allinds(index), setinds, compinds, opt.comps);
        [pcond pgroup pinter] = std_stat(erpdata, STUDY.etc.erpparams);
            
        if index == length(allinds), opt.legend = 'on'; end;
            std_plotcurve(alltimes, erpdata, 'condnames', STUDY.condition, 'legend', opt.legend, 'subject', opt.subject, ...
                                          'compinds', comp_names, 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottime, 'unitx', 'ms',  'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
        if length(allinds) > 1, 
            if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
            else                      title(sprintf('%s', opt.channels{index}));  
            end;
        end;
    end;
end;
