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
%   'plotmode'  - ['normal'|'condensed'|'none'] 'normal'  -> plot in a new figure; 
%                 'condensed' -> plot all curves in the current figure in a 
%                 condensed fashion. 'none' toggles off plotting {default: 'normal'}
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
% Known limitations: when plotting multiple clusters, the output
%                    contains the last plotted cluster.
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

function [STUDY, allersp, alltimes, allfreqs, pgroup, pcond, pinter] = std_erspplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erspstatplot;
    return;
end;

STUDY = pop_erspparams(STUDY, 'default');

[ opt moreparams ] = finputcheck( varargin, { ...
                               'design'      'integer' []              STUDY.currentdesign;
                               'topotime'    'real'    [] STUDY.etc.erspparams.topotime;
                               'topofreq'    'real'    [] STUDY.etc.erspparams.topofreq;
                               'maskdata'    'string'  [] 'off';
                               'timerange'   'real'    [] STUDY.etc.erspparams.timerange;
                               'freqrange'   'real'    [] STUDY.etc.erspparams.freqrange;
                               'ersplim'     'real'    [] STUDY.etc.erspparams.ersplim;
                               'caxis'       'real'    [] STUDY.etc.erspparams.ersplim;
                               'itclim'      'real'    [] STUDY.etc.erspparams.itclim;
                               'statistics'  'string'  [] STUDY.etc.erspparams.statistics;
                               'groupstats'  'string'  [] STUDY.etc.erspparams.groupstats;
                               'condstats'   'string'  [] STUDY.etc.erspparams.condstats;
                               'subbaseline' 'string'  [] STUDY.etc.erspparams.subbaseline;
                               'statmode'    'string'  [] ''; % deprecated
                               'threshold'   'real'    [] STUDY.etc.erspparams.threshold;
                               'naccu'       'integer' [] STUDY.etc.erspparams.naccu;
                               'mcorrect'    'string'  [] STUDY.etc.erspparams.mcorrect;
                               'singletrials' 'string' { 'on' 'off' }  STUDY.etc.erspparams.singletrials;
                               'channels'    'cell'    []              {};
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'itc' 'ersp' 'pac' } 'ersp';
                               'mode'        'string'  []              '';
                               'plottf'      'real'    []              [];
                               'comps'       {'integer','string'}  []              []; % for backward compatibility
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'plotmode'    'string' { 'normal' 'condensed' 'none' }  'normal';
                               'subject'     'string'  []              '' }, ...
                                  'std_erspstatplot', 'ignore');
if isstr(opt), error(opt); end;
if isempty(opt.caxis), 
    if strcmpi(opt.datatype, 'ersp')
         opt.caxis = opt.ersplim;
    elseif strcmpi(opt.datatype, 'itc')
        opt.caxis = [-opt.itclim opt.itclim];
    end;
end;

allconditions = STUDY.design(opt.design).variable(1).value;
allgroups     = STUDY.design(opt.design).variable(2).value;
paired = { fastif(strcmpi(STUDY.design(opt.design).variable(1).pairing, 'paired'), 'on', 'off') ...
           fastif(strcmpi(STUDY.design(opt.design).variable(2).pairing, 'paired'), 'on', 'off') };

% for backward compatibility
% --------------------------
if isempty(opt.plottf) & ~isempty(opt.topofreq) & ~isempty(opt.topotime) & ~isnan(opt.topofreq) & ~isnan(opt.topotime)
     opt.plottf = [ opt.topofreq(1) opt.topofreq(end) opt.topotime(1) opt.topotime(end) ];
end;
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if strcmpi(opt.singletrials, 'off') && ((~isempty(opt.subject) || ~isempty(opt.comps)))
    if strcmpi(opt.condstats, 'on') || strcmpi(opt.groupstats, 'on')
        opt.groupstats = 'off';
        opt.condstats   = 'off'; 
        disp('No statistics for single subject/component'); 
    end;
end;

if length(opt.comps) == 1
    opt.condstats = 'off'; opt.groupstats = 'off'; 
    disp('Statistics cannot be computed for single component');
end;

plottfopt = { ...
   'ersplim',     opt.caxis, ...
   'threshold',   opt.threshold, ...
   'maskdata',    opt.maskdata };
if ~isempty(opt.plottf) & length(opt.channels) < 5
    warndlg2(strvcat('ERSP/ITC parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;    

% plot single scalp map
% ---------------------
if ~isempty(opt.channels)

    [STUDY allersp alltimes allfreqs] = std_readersp(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', opt.datatype, 'subject', opt.subject, ...
        'singletrials', opt.singletrials, 'subbaseline', opt.subbaseline, 'timerange', opt.timerange, 'freqrange', opt.freqrange, 'design', opt.design);
    
    % select specific time and freq
    % -----------------------------
    if ~isempty(opt.plottf)
        if length(opt.plottf) < 3, 
            opt.plottf(3:4) = opt.plottf(2);
            opt.plottf(2) = opt.plottf(1);
        end;
        [tmp fi1] = min(abs(allfreqs-opt.plottf(1)));
        [tmp fi2] = min(abs(allfreqs-opt.plottf(2)));
        [tmp ti1] = min(abs(alltimes-opt.plottf(3)));
        [tmp ti2] = min(abs(alltimes-opt.plottf(4)));
        for index = 1:length(allersp(:))
            allersp{index} = mean(mean(allersp{index}(fi1:fi2,ti1:ti2,:,:),1),2);
            allersp{index} = reshape(allersp{index}, [1 size(allersp{index},3) size(allersp{index},4) ]);
        end;
        opt.plottf = { opt.plottf(1:2) opt.plottf(3:4) };
        [pcond pgroup pinter] = std_stat(allersp, 'groupstats', opt.groupstats, 'condstats', opt.condstats, 'paired', paired, ...
                                                  'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);
        if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end; % single subject STUDY                                
    else
        [pcond pgroup pinter] = std_stat(allersp, 'groupstats', opt.groupstats, 'condstats', opt.condstats, 'paired', paired, ...
                                                  'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);
        if (~isempty(pcond ) && (size( pcond{1},1) == 1 || size( pcond{1},2) == 1)) || ...
           (~isempty(pgroup) && (size(pgroup{1},1) == 1 || size(pgroup{1},2) == 1)), 
            pcond = {}; pgroup = {}; pinter = {}; 
            disp('No statistics possible for single subject STUDY');
        end; % single subject STUDY                                
    end
    
    % plot specific component
    % -----------------------
    if ~strcmpi(opt.plotmode, 'none')
        locs = eeg_mergelocs(ALLEEG.chanlocs);
        locs = locs(std_chaninds(STUDY, opt.channels));
        
        if ~isempty(opt.plottf)
            alltitles = std_figtitle('threshold', opt.threshold, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                     'statistics', opt.statistics, 'condnames', allconditions, 'cond2names', allgroups, 'chanlabels', { locs.labels }, ...
                                     'subject', opt.subject, 'valsunit', { 'Hz' 'ms' }, 'vals', opt.plottf, 'datatype', upper(opt.datatype));
            std_chantopo(allersp, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'caxis', opt.caxis, ...
                                          'chanlocs', locs, 'threshold', opt.threshold, 'titles', alltitles);
        else
            if length(opt.channels) > 1 & ~strcmpi(opt.plotmode, 'none'), figure; opt.plotmode = 'condensed'; end;
            nc = ceil(sqrt(length(opt.channels)));
            nr = ceil(length(opt.channels)/nc);
            for index = 1:max(cellfun(@(x)(size(x,3)), allersp(:)))
                if length(opt.channels) > 1, try, subplot(nr,nc,index, 'align'); catch, subplot(nr,nc,index); end; end;
                tmpersp = cell(size(allersp));
                for ind = 1:length(allersp(:))
                    if ~isempty(allersp{ind})
                        tmpersp{ind} = squeeze(allersp{ind}(:,:,index,:)); 
                    end;
                end;
                alltitles = std_figtitle('threshold', opt.threshold, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                         'statistics', opt.statistics, 'condnames', allconditions, 'cond2names', allgroups, 'chanlabels', { locs(index).labels }, ...
                                         'subject', opt.subject, 'datatype', upper(opt.datatype), 'plotmode', opt.plotmode);
                std_plottf(alltimes, allfreqs, tmpersp, 'datatype', opt.datatype, 'titles', alltitles, ...
                                           'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'plotmode', ...
                                           opt.plotmode, 'chanlocs', ALLEEG(1).chanlocs, plottfopt{:});
            end;
        end;
    end;
else
    
    if length(opt.clusters) > 1 & ~strcmpi(opt.plotmode, 'none'), figure; opt.plotmode = 'condensed'; end;
    nc = ceil(sqrt(length(opt.clusters)));
    nr = ceil(length(opt.clusters)/nc);
    comp_names = {};

    if length(opt.clusters) > 1 && ( strcmpi(opt.condstats, 'on') || strcmpi(opt.groupstats, 'on'))
        opt.condstats = 'off'; opt.groupstats = 'off';
    end;
    
    for index = 1:length(opt.clusters)

        [STUDY allersp alltimes allfreqs] = std_readersp(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'infotype', opt.datatype, ...
            'component', opt.comps, 'singletrials', opt.singletrials, 'subbaseline', opt.subbaseline, 'timerange', opt.timerange, 'freqrange', opt.freqrange, 'design', opt.design);
        if length(opt.clusters) > 1, try, subplot(nr,nc,index, 'align'); catch, subplot(nr,nc,index); end; end;

        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            comp_names = { STUDY.cluster(opt.clusters(index)).comps(opt.comps) };
            opt.subject = STUDY.datasetinfo(STUDY.cluster(opt.clusters(index)).sets(1,opt.comps)).subject;
        end;

        % select specific time and freq
        % -----------------------------
        if ~isempty(opt.plottf)
            if length(opt.plottf) < 3, 
                opt.plottf(3:4) = opt.plottf(2);
                opt.plottf(2) = opt.plottf(1);
            end;
            [tmp fi1] = min(abs(allfreqs-opt.plottf(1)));
            [tmp fi2] = min(abs(allfreqs-opt.plottf(2)));
            [tmp ti1] = min(abs(alltimes-opt.plottf(3)));
            [tmp ti2] = min(abs(alltimes-opt.plottf(4)));
            for index = 1:length(allersp(:))
                allersp{index} = mean(mean(allersp{index}(ti1:ti2,fi1:fi2,:,:),1),2);
                allersp{index} = reshape(allersp{index}, [1 size(allersp{index},3) size(allersp{index},4) ]);
            end;
            if opt.plottf(1) == opt.plottf(2), titlestr = [ num2str(opt.plottf(1)) ' Hz'];
            else                               titlestr = [ num2str(opt.plottf(1)) '-' num2str(opt.plottf(2)) ' Hz'];
            end;
            if opt.plottf(3) == opt.plottf(4), titlestr = [ ', ' num2str(opt.plottf(3)) ' ms'];
            else                               titlestr = [ ', ' num2str(opt.plottf(3)) '-' num2str(opt.plottf(4)) ' ms'];
            end;
            locs = eeg_mergelocs(ALLEEG.chanlocs);
            locs = locs(std_chaninds(STUDY, opt.channels));
        end

        [pcond pgroup pinter] = std_stat(allersp, 'groupstats', opt.groupstats, 'condstats', opt.condstats, 'paired', paired, ...
                                             'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);

        % plot specific component
        % -----------------------
        if index == length(opt.clusters), opt.legend = 'on'; end;
        if ~strcmpi(opt.plotmode, 'none')
            alltitles = std_figtitle('threshold', opt.threshold, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                     'statistics', opt.statistics, 'condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                     'subject', opt.subject, 'datatype', upper(opt.datatype), 'plotmode', opt.plotmode);
            
            std_plottf(alltimes, allfreqs, allersp, 'datatype', opt.datatype, ...
                                           'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'plotmode', ...
                                           opt.plotmode, 'titles', alltitles, ...
                                          'chanlocs', ALLEEG(1).chanlocs, plottfopt{:});
        end;
    end;
end;
