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

% $Log: not supported by cvs2svn $
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
if isempty(opt.plotfreq), opt.plotfreq = STUDY.etc.specparams.topofreq; end;
if isnan(  opt.plotfreq), opt.plotfreq = []; end;
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;

if ~isempty(opt.subject), groupstats = 'off'; disp('No group statistics for single subject');
else                      groupstats = STUDY.etc.specparams.groupstats;
end;
if ~isempty(opt.subject), condstats = 'off'; disp('No condition statistics for single subject');
else                      condstats = STUDY.etc.specparams.condstats;
end;
plotcurveopt = { ...
   'ylim',       STUDY.etc.specparams.ylim, ...
   'threshold',  STUDY.etc.specparams.threshold, ...
   'plotgroups',  STUDY.etc.specparams.plotgroups, ...
   'plotconditions', STUDY.etc.specparams.plotconditions, ...
   'statistics', STUDY.etc.specparams.statistics };

if ~isempty(opt.plotfreq) 
    if length(opt.channels) < 5
        warndlg2(strvcat('Spectrum parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
        return;
    end;
end;
if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'spec', 'freqrange', opt.freqrange, 'rmsubjmean', STUDY.etc.specparams.subtractsubjectmean);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'spec', 'freqrange', opt.freqrange, 'rmsubjmean', STUDY.etc.specparams.subtractsubjectmean);
end;

if strcmpi(opt.plotmode, 'condensed'), 
    plotcurveopt = { plotcurveopt{:} 'figure' 'off' }; 
end;
if length(allinds) > 1 & isempty(opt.channels)
    plotcurveopt = { plotcurveopt{:} 'figure' 'off' }; 
end;

% channel plotting
% ----------------
if ~isempty(opt.channels)
    structdat = STUDY.changrp;
    specdata = cell(size(structdat(allinds(1)).specdata));
    for ind =  1:length(structdat(allinds(1)).specdata(:))
        specdata{ind} = zeros([ size(structdat(allinds(1)).specdata{1}) length(allinds)]);
        for index = 1:length(allinds)
            specdata{ind}(:,:,index) = structdat(allinds(index)).specdata{ind};
            allfreqs                 = structdat(allinds(index)).specfreqs;
            compinds                 = structdat(allinds(index)).allinds;
            setinds                  = structdat(allinds(index)).setinds;
        end;
        specdata{ind} = squeeze(permute(specdata{ind}, [1 3 2])); % time elec subjects
    end;

    % select specific subject or component then plot
    % ----------------------------------------------
    if ~isempty(opt.subject), specdata = std_selsubject(specdata, opt.subject, setinds, { STUDY.datasetinfo(:).subject }, length(STUDY.subject)); end;
    
    % select specific time    
    % --------------------
    if ~isempty(opt.plotfreq)
        [tmp ti1] = min(abs(allfreqs-opt.plotfreq(1)));
        [tmp ti2] = min(abs(allfreqs-opt.plotfreq(end)));
        for index = 1:length(specdata(:))
            specdata{index} = mean(specdata{index}(ti1:ti2,:,:),1);
        end;
        if opt.plotfreq(1) == opt.plotfreq(end), titlestr = [ num2str(opt.plotfreq(1)) ' Hz'];
        else                                     titlestr = [ num2str(opt.plotfreq(1)) '-' num2str(opt.plotfreq(2)) ' ms'];
        end;
    end;
    
    % compute statistics and plot
    % ---------------------------
    [pcond pgroup pinter] = std_stat(specdata, STUDY.etc.specparams, 'groupstats', groupstats, 'condstats', condstats);
    locs = eeg_mergelocs(ALLEEG.chanlocs);
    locs = locs(std_chaninds(STUDY, opt.channels));
    if ~isempty(opt.plotfreq)
        std_chantopo(specdata, 'condnames', STUDY.condition, 'plottopo', fastif(length(allinds)==1, 'off', 'on'), ...
                                      'datatype', 'spec', 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'unitx', 'Hz', ...
                                      'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, 'topovals', titlestr, plotcurveopt{:});
    else
        std_plotcurve(allfreqs, specdata, 'condnames', STUDY.condition, 'plottopo', fastif(length(allinds)==1, 'off', 'on'), ...
                                      'datatype', 'spec', 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'unitx', 'Hz', ...
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

    for index = 1:length(allinds)

        if length(allinds) > 1, subplot(nr,nc,index); end;
        specdata  = STUDY.cluster(allinds(index)).specdata;
        allfreqs  = STUDY.cluster(allinds(index)).specfreqs;
        compinds  = STUDY.cluster(allinds(index)).allinds;
        setinds   = STUDY.cluster(allinds(index)).setinds;

        % plot specific component
        % -----------------------
        [specdata opt.subject comp_names] = std_selcomp(STUDY, specdata, allinds(index), setinds, compinds, opt.comps);
        [pcond pgroup pinter] = std_stat(specdata, STUDY.etc.specparams);
            
        if index == length(allinds), opt.legend = 'on'; end;
        std_plotcurve(allfreqs, specdata, 'condnames', STUDY.condition, 'legend', opt.legend, 'subject', opt.subject, ...
                                          'compinds', comp_names, 'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plotfreq, 'unitx', 'Hz',  'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, ...
                                          'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
        if length(allinds) > 1, 
            if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
            else                      title(sprintf('%s', opt.channels{index}));  
            end;
        end;
    end;
end;
