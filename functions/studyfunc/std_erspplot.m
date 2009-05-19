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
%   'plotmode'  - ['normal'|'condensed'] 'normal'  -> plot in a new figure; 
%                 'condensed' -> plot all curves in the current figure in a 
%                 condensed fashion {default: 'normal'}
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
% Revision 1.64  2009/05/13 00:52:12  arno
% fixing the ouput for scalp maps
%
% Revision 1.63  2008/04/16 17:56:33  arno
% fix color axis for scalp maps
%
% Revision 1.62  2007/11/15 04:01:36  arno
% *** empty log message ***
%
% Revision 1.61  2007/11/15 03:46:19  arno
% fix itclim
%
% Revision 1.60  2007/09/11 10:48:52  arno
% fix numerous display bugs
%
% Revision 1.59  2007/08/22 01:18:32  arno
% better handling of optional arguments
%
% Revision 1.58  2007/08/20 18:49:58  arno
% subtracting common baseline
%
% Revision 1.57  2007/08/14 02:51:30  arno
% fix statistics for several components
%
% Revision 1.56  2007/08/14 02:47:16  allen
% fix statistics for more than one component
%
% Revision 1.55  2007/08/14 02:09:20  arno
% allow plotting of individual components within the cluster
%
% Revision 1.54  2007/08/13 23:24:07  nima
% _
%
% Revision 1.53  2007/08/13 20:36:27  nima
% _
%
% Revision 1.52  2007/08/12 23:56:17  arno
% fix plotting topographic maps
%
% Revision 1.51  2007/08/12 23:33:53  arno
% fixing data ersp read if channel 1 is not included
%
% Revision 1.50  2007/08/10 22:42:53  arno
% no statistics for single components
%
% Revision 1.49  2007/08/10 19:47:56  nima
% _
%
% Revision 1.48  2007/06/25 07:19:38  toby
% altered multiple dipole plot windows to indicate number of components and subjects
%
% Revision 1.47  2007/05/11 02:12:16  toby
% statatistical mask wasn't plotting
%
% Revision 1.46  2007/05/01 21:17:04  arno
% clarify help message
%
% Revision 1.45  2007/04/30 20:57:27  arno
% update header and code for options
%
% Revision 1.44  2007/04/30 20:53:06  arno
% *** empty log message ***
%
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

[ opt moreparams ] = finputcheck( varargin, { ...
                               'topotime'    'real'    [] STUDY.etc.erspparams.topotime;
                               'topofreq'    'real'    [] STUDY.etc.erspparams.topofreq;
                               'maskdata'    'string'  [] STUDY.etc.erspparams.maskdata;
                               'timerange'   'real'    [] STUDY.etc.erspparams.timerange;
                               'freqrange'   'real'    [] STUDY.etc.erspparams.freqrange;
                               'ersplim'     'real'    [] STUDY.etc.erspparams.ersplim;
                               'itclim'      'real'    [] STUDY.etc.erspparams.itclim;
                               'statistics'  'string'  [] STUDY.etc.erspparams.statistics;
                               'groupstats'  'string'  [] STUDY.etc.erspparams.groupstats;
                               'condstats'   'string'  [] STUDY.etc.erspparams.condstats;
                               'subbaseline' 'string'  [] STUDY.etc.erspparams.subbaseline;
                               'statmode'    'string'  [] STUDY.etc.erspparams.statmode;
                               'threshold'   'real'    [] STUDY.etc.erspparams.threshold;
                               'naccu'       'integer' [] STUDY.etc.erspparams.naccu;
                               'channels'    'cell'    []              {};
                               'caxis'       'real'    []              [];
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'itc' 'ersp' } 'ersp';
                               'mode'        'string'  []              '';
                               'plottf'      'real'    []              [];
                               'comps'       {'integer','string'}  []              []; % for backward compatibility
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'subject'     'string'  []              '' }, ...
                                  'std_erspstatplot', 'ignore');
if isstr(opt), error(opt); end;
if isempty(opt.caxis), 
    if strcmpi(opt.datatype, 'ersp')
         opt.caxis = STUDY.etc.erspparams.ersplim;
    else opt.caxis = STUDY.etc.erspparams.itclim;
    end;
end;

% for backward compatibility
% --------------------------
if isempty(opt.plottf) & ~isempty(opt.topofreq) & ~isempty(opt.topotime)
     opt.plottf = [ opt.topofreq(1) opt.topofreq(end) opt.topotime(1) opt.topotime(end) ];
end;
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;

if ~isempty(opt.subject), opt.groupstats = 'off'; disp('No group statistics for single subject'); end;
if ~isempty(opt.subject), opt.condstats = 'off'; disp('No condition statistics for single subject'); end;

if length(opt.comps) == 1
    opt.condstats = 'off'; opt.groupstats = 'off'; 
    disp('Statistics cannot be computed for single component');
end;

plotcurveopt = { ...
   'ersplim',     fastif(strcmpi(opt.datatype, 'ITC'), [-opt.itclim opt.itclim], opt.ersplim), ...
   'threshold',   opt.threshold, ...
   'maskdata',    opt.maskdata, ...
   'groupstats',  opt.groupstats, ...
   'condstats',   opt.condstats, ...
   'statistics',  opt.statistics };
if ~isempty(opt.plottf) & length(opt.channels) < 5
    warndlg2(strvcat('ERSP/ITC parameters indicate that you wish to plot scalp maps', 'Select at least 5 channels to plot topography'));
    return;
end;    

% read data from disk
% -------------------
if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', opt.datatype, 'statmode', opt.statmode, 'subbaseline', opt.subbaseline);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', opt.datatype, 'statmode', opt.statmode, 'subbaseline', opt.subbaseline);
end;
opt.legend = 'off';

% plot single scalp map
% ---------------------
if ~isempty(opt.plottf)
    firstind = 1;
    while isempty(STUDY.changrp(firstind).erspdata)
        firstind = firstind+1;
    end;
    allersp = cell(size(STUDY.changrp(firstind).erspdata));
   
    for ind = 1:length(STUDY.changrp(firstind).erspdata(:))
        if size(STUDY.changrp(firstind).erspdata{1},3) == 1
             allersp{ind} = zeros([ size(STUDY.changrp(firstind).erspdata{1}) 1 length(opt.channels)]);
        else allersp{ind} = zeros([ size(STUDY.changrp(firstind).erspdata{1}) length(opt.channels)]);
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
                                      'chanlocs', locs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    % reselect time and frequency
    % ---------------------------
    maxind   = max(find( alltimes <= opt.plottf(4)));
    minind   = min(find( alltimes >= opt.plottf(3)));
    fmaxind  = max(find( allfreqs <= opt.plottf(2)));
    fminind  = min(find( allfreqs >= opt.plottf(1)));
    alltimes = mean(alltimes(minind:maxind));
    allfreqs = mean(allfreqs(fminind:fmaxind));
    for i =1:length(allersp(:))
        allersp{i}  = squeeze(mean(mean(allersp{i}(minind:maxind,fminind:fmaxind,:,:),1),2));
    end;
    return;
end;

if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);
comp_names = {};

for index = 1:length(allinds)

    if length(allinds) > 1, try, subplot(nr,nc,index, 'align'); catch, subplot(nr,nc,index); end; end;
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
            for ind = length(compinds{1,grpind}):-1:1
                if ~any(compinds{1,grpind}(ind) == comps) | ~any(setinds{1,grpind}(ind) == sets)
                    allersp{c}(:,:,ind) = [];
                else
                    comp_names{c,1} = comps;
                end;
            end;
        end;
        opt.subject = STUDY.datasetinfo(sets(1)).subject;
    end;
    
    % plot specific component
    % -----------------------
    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plottf(alltimes, allfreqs, allersp, 'condnames', STUDY.condition, 'subject', opt.subject, ...
                                       'legend', opt.legend, 'compinds', comp_names, 'datatype', opt.datatype,'plotmode', ...
                                       opt.plotmode, 'groupnames', STUDY.group, 'topovals', opt.plottf, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if isempty(opt.channels), %title(sprintf('Cluster %d', allinds(index)));
        if length(allinds) > 1, 
            title([ STUDY.cluster(allinds(index)).name ' (' num2str(length(STUDY.cluster(allinds(index)).comps)),' ICs, ' ...
                    num2str(length(unique(STUDY.cluster(allinds(index)).sets(1,:)))) ' Ss)' ]);            
            if strcmp(opt.datatype,'ersp')
                set(gcf,'name','Cluster ERSPs')
            elseif strcmp(opt.datatype,'itc')
                set(gcf,'name','Cluster ITCs')
            end;
        elseif ~strcmp(opt.mode,'together') % if it is not the mean ERSP that is being shown (which is the case when 'cluster properties' is plotted then put cluster number on the corner of figure
            h = gca;
            axes('position',[0.04 0.96 0.1 0.06]); 
            text(0,0,[STUDY.cluster(allinds(index)).name],'fontsize',13 );
            axis off;
            if strcmp(opt.datatype,'ersp')
                if length(opt.comps) ~= 1
                    set(gcf,'name',['ERSP of ' STUDY.cluster(allinds(index)).name])
                else
                    set(gcf,'name',['ERSP of a Component from cluster ' STUDY.cluster(allinds(index)).name])
                end;
            elseif strcmp(opt.datatype,'itc')
                if isempty(opt.comps)
                    set(gcf,'name',['ITC of ' STUDY.cluster(allinds(index)).name]);
                else
                    set(gcf,'name',['ITC of a Component from cluster ' STUDY.cluster(allinds(index)).name]);
                end;
            end;
            axes(h);
        end;
    else                      
        title(sprintf('%s', opt.channels{index}));  
    end;
end;
