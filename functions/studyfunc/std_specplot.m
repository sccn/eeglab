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
%  >> [STUDY] = std_specplot(STUDY, ALLEEG, key1, val1, key2, val2, ...);  
%  >> [STUDY specdata specfreqs pgroup pcond pinter] = std_specplot(STUDY, ALLEEG, ...);
%
% Inputs:
%   STUDY      - STUDY structure comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - vector of EEG dataset structures for the dataset(s) in the STUDY, 
%                typically created using load_ALLEEG().  
% Optional inputs for component plotting:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                spectrums of the requested clusters are plotted in the same figure, 
%                with spectrums for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, spectrum for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component spectrum plus the
%                average cluster spectrum in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel spectrum is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot spectrum of all subjects.
%
% Other optional inputs:
%   'plotmode'  - ['normal'|'condensed'] 'normal'  -> plot in a new figure; 
%                 'condensed' -> plot all curves in the current figure in a 
%                 condensed fashion {default: 'normal'}
%   'key','val' - All optional inputs to pop_specparams() are also accepted here
%                 to plot subset of time, statistics etc. The values used by default
%                 are the ones set using pop_specparams() and stored in the
%                 STUDY structure.
% Outputs:
%   STUDY      - the input STUDY set structure with the plotted cluster mean spectra
%                added?? to allow quick replotting.
%   specdata   - [cell] spectral data for each condition, group and subjects.
%                size of cell array is [nconds x ngroups]. Size of each element
%                is [freqs x subjects] for data channels or [freqs x components]
%                for component clusters. This array may be gicen as input 
%                directly to the statcond() function or std_stats() function
%                to compute statistics.
%   specfreqs  - [array] Sprectum point frequency values.
%   pgroup     - [array or cell] p-values group statistics. Output of the 
%                statcond() function.
%   pcond      - [array or cell] condition statistics. Output of the statcond() 
%                function.
%   pinter     - [array or cell] groups x conditions statistics. Output of
%                statcond() function.
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

function [STUDY, specdata, allfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_specplot;
    return;
end;

[STUDY, specdata, allfreqs, pgroup, pcond, pinter] = std_erpplot(STUDY, ALLEEG, 'datatype', 'spec', 'unitx', 'Hz', varargin{:});
return;

STUDY = pop_specparams(STUDY, 'default');

opt = finputcheck( varargin, { 'topofreq'    'real'    [] STUDY.etc.specparams.topofreq;
                               'freqrange'   'real'    [] STUDY.etc.specparams.freqrange;
                               'ylim'        'real'    [] STUDY.etc.specparams.ylim;
                               'caxis'       'real'    [] STUDY.etc.specparams.ylim;
                               'statistics'  'string'  [] STUDY.etc.specparams.statistics;
                               'groupstats'  'string'  [] STUDY.etc.specparams.groupstats;
                               'condstats'   'string'  [] STUDY.etc.specparams.condstats;
                               'plotgroups'  'string'  [] STUDY.etc.specparams.plotgroups;
                               'plotconditions'      'string' [] STUDY.etc.specparams.plotconditions;
                               'subtractsubjectmean' 'string' [] STUDY.etc.specparams.subtractsubjectmean
                               'threshold'   'real'    [] STUDY.etc.specparams.threshold;
                               'mcorrect'    'string'  [] STUDY.etc.specparams.mcorrect;
                               'naccu'       'integer' [] STUDY.etc.specparams.naccu;
                               'plotstderr'  'string'  []              'off';
                               'channels'    'cell'    []              {};
                               'clusters'    'integer' []              [];
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       { 'string' 'integer' } [] []; % for backward compatibility
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'singletrials' 'string' { 'on' 'off' }  'off';
                               'subject'     'string' []              '';
                               'statmode'    'string' { 'individual' 'common' 'trials' } 'individual'}, 'std_specplot');

if isstr(opt), error(opt); end;
if isstr(opt.comps), opt.comps = []; opt.plotsubjects = 'on'; end; % comps all

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if strcmpi(opt.condstats, 'on') || strcmpi(opt.groupstats, 'on') || ...
        (~isempty(opt.subject) || ~isempty(opt.comps)) && strcmpi(opt.singletrials, 'off'), 
    opt.groupstats = 'off';
    opt.constats   = 'off'; 
    disp('No statistics for single subject/component'); 
end;
plotcurveopt = { ...
   'ylim',           opt.ylim, ...
   'threshold',      opt.threshold, ...
   'plotgroups',     opt.plotgroups, ...
   'plotconditions', opt.plotconditions  };

if ~isnan(opt.topofreq) & length(opt.channels) < 5
    warndlg2(strvcat('Spectrum parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
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
    [STUDY specdata allfreqs setinds allinds] = std_readspec(STUDY, ALLEEG, 'channels', opt.channels, 'freqrange', opt.freqrange, ...
        'rmsubjmean', opt.subtractsubjectmean, 'subject', opt.subject, 'singletrials', opt.singletrials);
    if isempty(specdata), return; end;
        
    % select specific time    
    % --------------------
    if ~isempty(opt.topofreq) & ~isnan(opt.topofreq)
        [tmp ti1] = min(abs(allfreqs-opt.topofreq(1)));
        [tmp ti2] = min(abs(allfreqs-opt.topofreq(end)));
        for index = 1:length(specdata(:))
            specdata{index} = mean(specdata{index}(ti1:ti2,:,:),1);
        end;
    end;
    
    % compute statistics and plot
    % ---------------------------
    [pcond pgroup pinter] = std_stat(specdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, ...
                                         'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect );
    locs = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locs(std_chaninds(STUDY, opt.channels));
    alltitles = std_figtitle('threshold', opt.threshold, 'mcorrect', opt.mcorrect, 'plotsubjects', opt.plotsubjects, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                             'statistics', opt.statistics, 'condnames', STUDY.condition, 'cond2names', STUDY.group, 'chanlabels', { locs.labels }, ...
                             'subject', opt.subject, 'valsunit', 'Hz', 'vals', opt.topofreq, 'datatype', 'Spectrum', 'cond2group', opt.plotgroups, 'condgroup', opt.plotconditions);
    
    if ~isempty(opt.topofreq) & all(~isnan(opt.topofreq))
        std_chantopo(specdata, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'caxis', opt.caxis, ...
                                      'chanlocs', locs, 'threshold', opt.threshold, 'titles', alltitles);
    else
        std_plotcurve(allfreqs, specdata, 'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
            'chanlocs', locs, 'titles', alltitles, 'plotsubjects', opt.plotsubjects, 'unitx', 'Hz', 'plotstderr', opt.plotstderr, ...
            'condnames', STUDY.condition, 'groupnames', STUDY.group, 'plottopo', fastif(length(allinds) > 5, 'on', 'off'), plotcurveopt{:});
    end;
    set(gcf,'name','Channel Spectra');
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
        
        [STUDY specdata allfreqs setinds compinds] = std_readspec(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'freqrange', opt.freqrange, ...
                    'rmsubjmean', opt.subtractsubjectmean, 'component', opt.comps, 'singletrials', opt.singletrials);
        if isempty(specdata), return; end;
        
        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            comp_names = { STUDY.cluster(opt.clusters(index)).comps(opt.comps) };
            opt.subject = STUDY.datasetinfo(STUDY.cluster(opt.clusters(index)).sets(1,opt.comps)).subject;
        end;
        [pcond pgroup pinter] = std_stat(specdata, 'groupstats', opt.groupstats, 'condstats', opt.condstats, ...
                                         'statistics', opt.statistics, 'naccu', opt.naccu, 'threshold', opt.threshold, 'mcorrect', opt.mcorrect);
            
        if index == length(opt.clusters), opt.legend = 'on'; end;
        alltitles = std_figtitle('threshold', opt.threshold, 'plotsubjects', opt.plotsubjects, 'mcorrect', opt.mcorrect, 'condstat', opt.condstats, 'cond2stat', opt.groupstats, ...
                                 'statistics', opt.statistics, 'condnames', STUDY.condition, 'cond2names', STUDY.group, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                 'subject', opt.subject, 'datatype', 'spectrum', 'cond2group', opt.plotgroups, 'condgroup', opt.plotconditions);
        
        std_plotcurve(allfreqs, specdata, 'condnames', STUDY.condition, 'legend', opt.legend, 'groupnames', STUDY.group, 'plotstderr', opt.plotstderr, ...
                                          'titles', alltitles, 'unitx', 'Hz',  'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    end;
    set(gcf,'name','Component Spectra');
    axcopy(gca);
end;
