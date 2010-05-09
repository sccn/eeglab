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

% $Log: std_specplot.m,v $
% Revision 1.62  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%
% Revision 1.61  2010/02/09 06:07:27  arno
% Fixed new title problem and implemented 3-level significance
%
% Revision 1.60  2010/02/06 05:47:53  arno
% New titles for figures
%
% Revision 1.59  2009/10/07 05:07:19  arno
% Fix missing conditions/groups
%
% Revision 1.58  2009/08/29 04:24:56  arno
% new statistics
%
% Revision 1.56  2009/08/29 00:38:32  arno
% move all statistics to std_stat
%
% Revision 1.55  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.54  2008/09/25 15:04:18  arno
% allow plotting groups of different sizes
%
% Revision 1.53  2008/04/16 17:55:27  arno
% fix color axis for scalp maps
%
% Revision 1.51  2007/09/11 10:55:03  arno
% fix numerous small display bugs
%
% Revision 1.50  2007/08/25 01:05:01  arno
% nothing
%
% Revision 1.49  2007/08/14 19:20:26  nima
% _
%
% Revision 1.48  2007/08/13 23:24:20  nima
% _
%
% Revision 1.47  2007/08/10 22:33:56  arno
% no statistics for single components
%
% Revision 1.46  2007/08/09 18:32:58  arno
% only display warning message if statistics are set
%
% Revision 1.45  2007/08/09 18:29:59  arno
% add warning message and disable statistics for multiple cluster plotting
%
% Revision 1.44  2007/06/25 04:38:56  toby
% altered multiple dipole plot windows to indicate number of components and subjects
%
% Revision 1.43  2007/05/01 21:17:05  arno
% clarify help message
%
% Revision 1.42  2007/04/30 20:52:43  arno
% header
%
% Revision 1.41  2007/04/30 20:11:27  arno
% update header and code
%
% Revision 1.40  2007/04/28 00:28:11  arno
% fix backward compatibility
% .,
%
% Revision 1.39  2007/04/06 18:52:43  arno
% figure creation
%
% Revision 1.38  2007/04/05 22:02:05  arno
% plot in condensed mode
%
% Revision 1.37  2007/03/14 01:15:39  arno
% plot condensed mode
%
% Revision 1.36  2007/03/14 01:01:13  arno
% plotting spectrum of several clusters
%
% Revision 1.35  2007/02/28 12:03:58  arno
% output statistics and documentation
%
% Revision 1.34  2007/01/26 18:04:55  arno
% reprogrammed from scratch again
%
% Revision 1.33  2006/11/23 00:50:01  arno
% rmsubjmean flag
%
% Revision 1.32  2006/11/23 00:27:40  arno
% add subject info for components
%
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
            'chanlocs', locs, 'titles', alltitles, 'plotsubjects', opt.plotsubjects, 'unitx', 'Hz', ...
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
        
        std_plotcurve(allfreqs, specdata, 'condnames', STUDY.condition, 'legend', opt.legend, 'groupnames', STUDY.group,  ...
                                          'titles', alltitles, 'unitx', 'Hz',  'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    end;
    set(gcf,'name','Component Spectra');
    axcopy(gca);
end;
