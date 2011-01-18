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

allconditions = STUDY.design(opt.design).variable(1).value;
allgroups     = STUDY.design(opt.design).variable(2).value;
paired = { STUDY.design(opt.design).variable(1).pairing ...
           STUDY.design(opt.design).variable(2).pairing };

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

if ~isnan(opt.topotime) & length(opt.channels) < 5
    warndlg2(strvcat('ERP parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;

plotcurveopt = {};
if length(opt.clusters) > 1 || length(opt.channels) > 1
    plotcurveopt = { 'figure' 'off' }; 
    opt.plotconditions = 'together';
    opt.plotgroups     = 'together';
    opt.condstats  = 'off'; 
    opt.groupstats = 'off';
end;
plotcurveopt = { plotcurveopt{:} ...
   'ylim',           opt.ylim, ...
   'threshold',      opt.threshold, ...
   'unitx'           opt.unitx, ...
   'filter',         opt.filter, ...
   'plotgroups',     opt.plotgroups, ...
   'plotconditions', opt.plotconditions };

% channel plotting
% ----------------
if ~isempty(opt.channels)
    % plot channel
    % --------------
    if length(opt.channels) > 1, figure('color', 'w'); end;
    nc = ceil(sqrt(length(opt.channels)));
    nr = ceil(length(opt.channels)/nc);
    comp_names = {};
    
    % channel indices (for topomap, have to be transposed)
    chaninds = 1:length(opt.channels);
    if ~isempty(opt.topotime) && all(~isnan(opt.topotime))
        chaninds = chaninds';
    end;
    
    for index = chaninds
        
        if length(opt.channels) > 1, subplot(nr,nc,index); end;
        if strcmpi(opt.datatype, 'erp')
            [STUDY erpdata alltimes] = std_readerp(STUDY, ALLEEG, 'channels', opt.channels(index), 'timerange', opt.timerange, ...
                    'subject', opt.subject, 'singletrials', opt.singletrials, 'design', opt.design);
        else
            [STUDY erpdata alltimes] = std_readspec(STUDY, ALLEEG, 'channels', opt.channels(index), 'freqrange', opt.freqrange, ...
                'rmsubjmean', opt.subtractsubjectmean, 'subject', opt.subject, 'singletrials', opt.singletrials, 'design', opt.design);
        end;
        if isempty(erpdata), return; end;

        % select specific time    
        % --------------------
        if ~isempty(opt.topotime) & ~isnan(opt.topotime)
            [tmp ti1] = min(abs(alltimes-opt.topotime(1)));
            [tmp ti2] = min(abs(alltimes-opt.topotime(end)));
            for condind = 1:length(erpdata(:))
                if ~isempty(erpdata{condind})
                    erpdata{condind} = mean(erpdata{condind}(ti1:ti2,:,:),1);
                end;
            end;
        end;

        % compute statistics
        % ------------------
        [pcond pgroup pinter] = std_stat(erpdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, 'paired', paired, ...
                                             'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);
        if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end; % single subject STUDY                                
        if length(opt.channels) > 5 && ndims(erpdata{1}) < 3, pcond = {}; pgroup = {}; pinter = {}; end; % topo plotting for single subject

        % plot
        % ----
        locs = eeg_mergelocs(ALLEEG.chanlocs);
        locs = locs(std_chaninds(STUDY, opt.channels(index)));
        [alltitles alllegends ] = std_figtitle('threshold', opt.threshold, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                 'statistics', opt.statistics, 'condnames', allconditions, 'plotsubjects', opt.plotsubjects, 'cond2names', allgroups, 'chanlabels', { locs.labels }, ...
                                 'subject', opt.subject, 'valsunit', opt.unitx, 'vals', opt.topotime, 'datatype', datatypestr, 'cond2group', opt.plotgroups, 'condgroup', opt.plotconditions);

        if ~isempty(opt.topotime) && all(~isnan(opt.topotime))
            std_chantopo(erpdata, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'caxis', opt.caxis, ...
                                          'chanlocs', locs, 'threshold', opt.threshold, 'titles', alltitles);
        else
            if index < length(opt.channels), alllegends = {}; end;
            std_plotcurve(alltimes, erpdata, 'groupstats', pgroup, 'legend', alllegends, 'condstats', pcond, 'interstats', pinter, ...
                'chanlocs', locs, 'titles', alltitles, 'plotsubjects', opt.plotsubjects, ...
                'condnames', allconditions, 'groupnames', allgroups, 'plottopo', fastif(length(opt.channels) > 5, 'on', 'off'), plotcurveopt{:});
        end;
    end;
    
    set(gcf,'name',['Channel ' datatypestr ]);
    axcopy(gca);
else 
    % plot component
    % --------------
    if length(opt.clusters) > 1, figure('color', 'w'); end;
    nc = ceil(sqrt(length(opt.clusters)));
    nr = ceil(length(opt.clusters)/nc);
    comp_names = {};
    
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
            
        [alltitles alllegends ] = std_figtitle('threshold', opt.threshold, 'plotsubjects', opt.plotsubjects, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                 'statistics', opt.statistics, 'condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                 'subject', opt.subject, 'valsunit', opt.unitx, 'vals', opt.topotime, 'datatype', datatypestr, 'cond2group', opt.plotgroups, 'condgroup', opt.plotconditions);
        
        if index == length(opt.clusters), alllegends = {}; end;
        std_plotcurve(alltimes, erpdata, 'condnames', allconditions, 'legend', alllegends, 'groupnames', allgroups, ...
                                          'titles', alltitles, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    end;
    
    set(gcf,'name', ['Component ' datatypestr ] );
    axcopy(gca);
end;

