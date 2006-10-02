% std_erpplot() - Commandline function to plot STUDY cluster component ERPs. Either 
%                 displays mean ERP of all requested clusters in the same figure, with 
%                 ERPs for different conditions (if any) plotted in different colors. 
%                 Else, displays ERP for each specified cluster in separate figures 
%                 (per condition), each containing the cluster component ERPs plus 
%                 the grand mean cluster ERP (in bold). ERPs can be plotted only if 
%                 component ERPs were computed and saved in the STUDY EEG datasets. 
%                 These can be computed during pre-clustering using the gui-based 
%                 function pop_preclust() or the equivalent commandline functions 
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
%   'clusters' - [numeric vector]  -> specific cluster indices to plot.
%                            'all' -> plot all clusters in STUDY {default: 'all'}.
%                If no component index (next option) is given, the average ERPs of 
%                the requested clusters are plotted in the same figure, with ERPs for
%                different conditions (if any) plotted in different colors. In 'comps' 
%                mode, ERPS for each specified cluster are plotted in separate figures 
%                (per condition), each containing cluster component ERPs plus the
%                average cluster ERP in bold. Note this parameter has no effect 
%                if the 'comps' option is used. {default: 'centroid'}.
%   'comps'    - [numeric vector]  -> indices of the cluster components to plot.
%                            'all' -> plot all the components in the cluster 
%                {default: 'all'}.
%
% Optional inputs for channel plotting:
%   'changrp'  - [numeric vector]  -> specific channel group to plot. By
%                default, the channel grand ERP is plotted (using the same
%                format as the cluster component centroid described above)
%   'subject'  - [numeric vector]  -> index of the subject to plot for changrp.
%                            'all' -> plot all the components in the cluster 
%
% Other optional inputs:
%   'figure'   - ['on'|'off'] for the 'centroid' mode option. 
%                 'on'  -> plot in a new figure; 
%                 'off' -> plot in the current figure {default: 'on'}
% Outputs:
%   STUDY      - the input STUDY set structure modified with plotted cluster 
%                 mean ERP to allow quick replotting (unless cluster means 
%                 already exists in the STUDY).  
%   Example:
%              >> [STUDY] = std_erpplot(STUDY,ALLEEG, 'clusters', 2, 'comps', 'all');
%                 % Plot cluster-2 components ERPs plus the mean ERP in bold. 
%
%  See also  pop_clustedit(), pop_preclust()
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
                               'statmode'    'string'  { 'individual' 'common' 'trials' } 'individual'}, 'std_erpplot');
if isstr(opt), error(opt); end;

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if ~isempty(opt.comps), 
    opt.subject = STUDY.datasetinfo( STUDY.cluster(opt.clusters).sets(1,opt.comps)).subject;
end;

if ~isempty(opt.subject), statgroup = 'off'; disp('No group statistics for single subject');
else                      statgroup = STUDY.etc.erpparams.statgroup;
end;
if ~isempty(opt.subject), statcond = 'off'; disp('No condition statistics for single subject');
else                      statcond = STUDY.etc.erpparams.statcond;
end;
plotcurveopt = { ...
   'ylim',       STUDY.etc.erpparams.ylim, ...
   'threshold',  STUDY.etc.erpparams.threshold, ...
   'statgroup',  statgroup, ...
   'statcond',   statcond, ...
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
    [pgroup pcond pinter] = std_plot(alltimes, erpdata, 'condname', STUDY.condition, 'legend', opt.legend, ...
                                      'plotmode', opt.plotmode, 'groupname', STUDY.group, 'plotx', opt.plottime, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
