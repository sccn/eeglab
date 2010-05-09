% std_erpplot() - Command line function to plot STUDY cluster component ERPs. Either 
%                 displays grand mean ERPs for all requested clusters in the same figure, 
%                 with ERPs for different conditions (if any) plotted in different colors. 
%                 Else, displays ERP for each specified cluster in separate figures 
%                 (per condition), each containing the cluster component ERPs plus 
%                 the grand mean cluster ERP (in bold). ERPs can be plotted only if 
%                 component ERPs were computed and saved in the STUDY EEG
%                 datasets. 
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

% $Log: std_erpplot.m,v $
% Revision 1.68  2010/02/26 10:58:31  claire
% fixing unitx and disabling condstats and groupstats for single subjects
%
% Revision 1.67  2010/02/24 10:52:36  arno
% Implemented new single trial statistics
%
% Revision 1.66  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%
% Revision 1.65  2010/02/09 06:07:27  arno
% Fixed new title problem and implemented 3-level significance
%
% Revision 1.64  2010/02/06 05:47:52  arno
% New titles for figures
%
% Revision 1.63  2009/10/07 05:07:19  arno
% Fix missing conditions/groups
%
% Revision 1.62  2009/08/29 04:24:56  arno
% new statistics
%
% Revision 1.61  2009/07/02 19:05:08  arno
% plotting groups of unequal size
%
% Revision 1.60  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.59  2008/04/16 17:55:53  arno
% fix color axis for scalp maps
%
% Revision 1.57  2007/08/14 19:29:47  nima
% _
%
% Revision 1.56  2007/08/13 23:25:16  nima
% _
%
% Revision 1.55  2007/08/10 22:43:07  arno
% no statistics for single components
%
% Revision 1.54  2007/08/09 18:33:00  arno
% only display warning message if statistics are set
%
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
erpdata = []; alltimes = []; 
pgroup = []; pcond = []; pinter = [];

% find datatype and default options
% ---------------------------------
dtype = 'erp';
for ind = 1:2:length(varargin)
    if strcmpi(varargin{ind}, 'datatype')
        dtype = varargin{ind+1}; 
    end;
end;
eval( [ 'STUDY = pop_' dtype 'params(STUDY, ''default'');' ...
        'params = STUDY.etc.' dtype 'params;' ] );
fields     = { 'filter' 'subtractsubjectmean' 'timerange' 'freqrange' 'topotime' 'topofreq' };
defaultval = { [] 'off' [] [] [] [] };
for ind=1:length(fields)
    if ~isfield(params, fields{ind}), 
        params = setfield(params, fields{ind}, defaultval{ind}); 
    end;
end;

opt = finputcheck( varargin, { 'topotime'    'real'    [] params.topotime;
                               'topofreq'    'real'    [] params.topofreq;
                               'filter'      'real'    [] params.filter;
                               'timerange'   'real'    [] params.timerange;
                               'freqrange'   'real'    [] params.freqrange;
                               'ylim'        'real'    [] params.ylim;
                               'caxis'       'real'    [] params.ylim;                               
                               'statistics'  'string'  [] params.statistics;
                               'groupstats'  'string'  [] params.groupstats;
                               'condstats'   'string'  [] params.condstats;
                               'mcorrect'    'string'  [] params.mcorrect;
                               'plotgroups'  'string'  [] params.plotgroups;
                               'plotconditions' 'string'  [] params.plotconditions;
                               'subtractsubjectmean' 'string' [] params.subtractsubjectmean;
                               'threshold'   'real'    [] params.threshold;
                               'naccu'       'integer' [] params.naccu;
                               'singletrials' 'string' { 'on' 'off' }  params.singletrials;
                               'design'      'integer' []              STUDY.currentdesign;
                               'channels'    'cell'    []              {};
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'erp' 'spec' } 'erp';      
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       { 'string' 'integer' } [] []; % for backward compatibility
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'unitx'       'string' { 'ms' 'Hz' }    'ms';
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'subject'     'string'  []              '';
                               'statmode'    'string'  { 'subjects' 'common' 'trials' } 'subjects'}, 'std_erpplot');
if isstr(opt), error(opt); end;
if isstr(opt.comps), opt.comps = []; opt.plotsubjects = 'on'; end;
if ~isempty(opt.topofreq),  opt.topotime  = opt.topofreq; end;
if ~isempty(opt.freqrange), opt.timerange = opt.freqrange; end;
datatypestr = upper(opt.datatype);
if strcmpi(datatypestr, 'spec'), datatypestr = 'Spectrum'; end;

% =======================================================================
% below this line, all the code should be non-specific to ERP or spectrum
% =======================================================================

allconditions = STUDY.design(opt.design).condition;
allgroups     = STUDY.design(opt.design).group;
paired = { fastif(strcmpi(STUDY.design(opt.design).statvar1, 'paired'), 'on', 'off') ...
           fastif(strcmpi(STUDY.design(opt.design).statvar2, 'paired'), 'on', 'off') };

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if strcmpi(opt.singletrials, 'off') && ((~isempty(opt.subject) || ~isempty(opt.comps)))
    if strcmpi(opt.condstats, 'on') || strcmpi(opt.groupstats, 'on')
        opt.groupstats = 'off';
        opt.condstats   = 'off'; 
        disp('No statistics for single subject/component'); 
    end;
end;
plotcurveopt = { ...
   'ylim',           opt.ylim, ...
   'threshold',      opt.threshold, ...
   'unitx'           opt.unitx, ...
   'filter',         opt.filter, ...
   'plotgroups',     opt.plotgroups, ...
   'plotconditions', opt.plotconditions };

if ~isnan(opt.topotime) & length(opt.channels) < 5
    warndlg2(strvcat('ERP parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;

if length(opt.clusters) > 1 && isempty(opt.channels)
    plotcurveopt = { plotcurveopt{:} 'figure' 'off' }; 
    opt.plotconditions = 'together';
    opt.plotgroups     = 'together';
end;

% channel plotting
% ----------------
if ~isempty(opt.channels)
    if strcmpi(opt.datatype, 'erp')
        [STUDY erpdata alltimes] = std_readerp(STUDY, ALLEEG, 'channels', opt.channels, 'timerange', opt.timerange, ...
                'subject', opt.subject, 'singletrials', opt.singletrials, 'design', opt.design);
    else
        [STUDY erpdata alltimes] = std_readspec(STUDY, ALLEEG, 'channels', opt.channels, 'freqrange', opt.freqrange, ...
            'rmsubjmean', opt.subtractsubjectmean, 'subject', opt.subject, 'singletrials', opt.singletrials, 'design', opt.design);
    end;
    if isempty(erpdata), return; end;
    
    % select specific time    
    % --------------------
    if ~isempty(opt.topotime) & ~isnan(opt.topotime)
        [tmp ti1] = min(abs(alltimes-opt.topotime(1)));
        [tmp ti2] = min(abs(alltimes-opt.topotime(end)));
        for index = 1:length(erpdata(:))
            erpdata{index} = mean(erpdata{index}(ti1:ti2,:,:),1);
        end;
    end;
    
    % compute statistics and plot
    % ---------------------------
    [pcond pgroup pinter] = std_stat(erpdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, 'paired', paired, ...
                                         'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);
    if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end; % single subject STUDY                                
    locs = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locs(std_chaninds(STUDY, opt.channels));
    alltitles = std_figtitle('threshold', opt.threshold, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                             'statistics', opt.statistics, 'condnames', allconditions, 'plotsubjects', opt.plotsubjects, 'cond2names', allgroups, 'chanlabels', { locs.labels }, ...
                             'subject', opt.subject, 'valsunit', opt.unitx, 'vals', opt.topotime, 'datatype', datatypestr, 'cond2group', opt.plotgroups, 'condgroup', opt.plotconditions);
    
    if ~isempty(opt.topotime) && all(~isnan(opt.topotime))
        std_chantopo(erpdata, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'caxis', opt.caxis, ...
                                      'chanlocs', locs, 'threshold', opt.threshold, 'titles', alltitles);
    else
        std_plotcurve(alltimes, erpdata, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
            'chanlocs', locs, 'titles', alltitles, 'plotsubjects', opt.plotsubjects, ...
            'condnames', allconditions, 'groupnames', allgroups, 'plottopo', fastif(length(opt.channels) > 5, 'on', 'off'), plotcurveopt{:});
    end;

    set(gcf,'name',['Channel ' datatypestr ]);
else 
    % plot component
    % --------------
    opt.legend = 'off';
    if length(opt.clusters) > 1, figure('color', 'w'); end;
    nc = ceil(sqrt(length(opt.clusters)));
    nr = ceil(length(opt.clusters)/nc);
    comp_names = {};

    if length(opt.clusters) > 1 && ( strcmpi(opt.condstats, 'on') || strcmpi(opt.groupstats, 'on'))
        opt.condstats = 'off'; opt.groupstats = 'off';
    end;
    
    for index = 1:length(opt.clusters)

        if length(opt.clusters) > 1, subplot(nr,nc,index); end;
        if strcmpi(opt.datatype, 'erp')
            [STUDY erpdata alltimes] = std_readerp(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'timerange', opt.timerange, ...
                        'component', opt.comps, 'singletrials', opt.singletrials, 'design', opt.design);
        else
            [STUDY erpdata alltimes] = std_readspec(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'freqrange', opt.freqrange, ...
                        'rmsubjmean', opt.subtractsubjectmean, 'component', opt.comps, 'singletrials', opt.singletrials, 'design', opt.design);
        end;
        if isempty(erpdata), return; end;

        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            comp_names = { STUDY.cluster(opt.clusters(index)).comps(opt.comps) };
            opt.subject = STUDY.datasetinfo(STUDY.cluster(opt.clusters(index)).sets(1,opt.comps)).subject;
        end;
        [pcond pgroup pinter] = std_stat(erpdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, 'paired', paired, ...
                                         'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);
            
        if index == length(opt.clusters), opt.legend = 'on'; end;
        alltitles = std_figtitle('threshold', opt.threshold, 'plotsubjects', opt.plotsubjects, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                 'statistics', opt.statistics, 'condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                 'subject', opt.subject, 'valsunit', opt.unitx, 'vals', opt.topotime, 'datatype', datatypestr, 'cond2group', opt.plotgroups, 'condgroup', opt.plotconditions);
        
        std_plotcurve(alltimes, erpdata, 'condnames', allconditions, 'legend', opt.legend, 'groupnames', allgroups,  ...
                                          'titles', alltitles, 'unitx', opt.unitx,  'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    end;
    
    set(gcf,'name', ['Component ' datatypestr ] );
    axcopy(gca);
end;

