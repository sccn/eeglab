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
%              >> [STUDY] = std_erpplot(STUDY, ALLEEG, key1, val1, key2, val2);  
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
%   'changrp'  - [numeric vector]  specific channel group to plot. By
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
%   Example:
%            >> [STUDY] = std_erpplot(STUDY,ALLEEG, 'clusters', 2, 'comps', 'all');
%               % Plot cluster-2 component ERPs plus the mean ERP in bold. 
%
%  See also  pop_clustedit(), pop_preclust(), eeg_createdata(), eeg_preclust(). std_propplot()
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
                            
function [STUDY allerp alltimes ] = std_erpplot(STUDY, ALLEEG, varargin)

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
if ~isempty(opt.comps), 
    opt.subject = STUDY.datasetinfo( STUDY.cluster(opt.clusters).sets(1,opt.comps)).subject;
end;

if ~isempty(opt.subject), groupstats = 'off'; disp('No group statistics for single subject');
else                      groupstats = STUDY.etc.erpparams.groupstats;
end;
if ~isempty(opt.subject), condstats = 'off'; disp('No condition statistics for single subject');
else                      condstats = STUDY.etc.erpparams.condstats;
end;
plotcurveopt = { ...
   'ylim',       STUDY.etc.erpparams.ylim, ...
   'threshold',  STUDY.etc.erpparams.threshold, ...
   'groupstats',  groupstats, ...
   'condstats',   condstats, ...
   'plotgroups',  STUDY.etc.erpparams.plotgroups, ...
   'plotconditions',   STUDY.etc.erpparams.plotconditions, ...
   'statistics', STUDY.etc.erpparams.statistics };

if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'erp', 'timerange', opt.timerange);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'erp', 'timerange', opt.timerange);
end;

opt.legend = 'off';
if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);

for index = 1:length(allinds)

    if length(allinds) > 1, subplot(nr,nc,index); end;
    if ~isempty(opt.channels)
        erpdata  = STUDY.changrp(allinds(index)).erpdata;
        alltimes = STUDY.changrp(allinds(index)).erptimes;
    else
        erpdata  = STUDY.cluster(allinds(index)).erpdata;
        alltimes = STUDY.cluster(allinds(index)).erptimes;
    end;

    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(erpdata,1)
            for g = 1:size(erpdata,2)
                erpdata{c,g} = erpdata{c,g}(:,subjind);
            end;
        end;
    end;

    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plot(alltimes, erpdata, 'condnames', STUDY.condition, 'legend', opt.legend, ...
                                      'plotmode', opt.plotmode, 'groupnames', STUDY.group, 'plottopo', opt.plottime, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
