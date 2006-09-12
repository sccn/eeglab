% std_erspplot() - plot cluster ERSPs. Displays either mean cluster ERSPs, 
%                  or else all cluster component ERSPs plus the mean cluster 
%                  ERSP in one figure per condition. The ERSPs can be plotted 
%                  only if component ERSPs were computed and saved in the 
%                  EEG datasets in the STUDY. These may either be computed 
%                  during pre-clustering using the gui-based function 
%                  pop_preclust(), or via the equivalent commandline functions 
%                  eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%        >> [STUDY] = std_erspplot(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [numeric vector]  -> cluster numbers to plot.
%                            'all' -> plot all clusters in STUDY 
%                            {default: 'all'}.
%   'comps'    - [numeric vector]  -> cluster components to plot.
%                            'all' -> plot all cluster components 
%                            {default: 'all'}.
%   'channels' - [numeric vector]  -> channels to plot.
%   'mode'     - ['centroid'|'individual'] plotting mode. In 'centroid' 
%                mode, the average ERSPs of the requested clusters or channels  
%                are plotted in the same figure - one per condition. In 
%                'individual' mode, component ERSPs for each
%                cluster (or channel) are plotted in a separate 
%                figure (per condition) with the mean ERSP. 
%                Note that for clusters, this option is irrelevant if component  
%                indices are provided as input {default: 'centroid'} 
%   'figure'  - ['on'|'off'] 'on' -> plot on a new figure; 
%                'off' -> plot on current figure 'figure'.
%                Note: 'off' is optional for one cluster in 'centroid' mode.
%                Useful for incomporating a cluster ERSP into 
%                a complex figure. In case of multiple conditions, 
%                only the first condition is displayed, but clicking on 
%                the figure will open a new figure with all conditions 
%                plotted separately {default: 'on'} 
% Output:
%   STUDY     - the input STUDY set structure modified with plotted cluster 
%               mean ERSPs to allow quick replotting (unless cluster means 
%               already exists in the STUDY).  
%
% Example:
%        >> [STUDY] = std_erspplot(STUDY,ALLEEG, 'clusters', 'all', ...
%                                       'mode', 'centroid');
%           % Plot the mean ERSPs of all clusters in STUDY on the same figure. 
%
% See also: pop_clustedit(), pop_preclust()
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
                            
function [STUDY allersp alltimes ] = std_erspstatplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erspstatplot;
    return;
end;

STUDY = pop_erspparams(STUDY, 'default');

opt = finputcheck( varargin, { 'channels'    'cell'    []              {};
                               'caxis'       'real'    []              [];
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'itc' 'ersp' } 'ersp';
                               'mode'        'string'  []              '';
                               'plottf'      'real'    []              [];
                               'timerange'   'real'    []              STUDY.etc.erspparams.timerange;
                               'freqrange'   'real'    []              STUDY.etc.erspparams.freqrange;
                               'mode'        'string'  []              ''; % for backward compatibility
                               'comps'       'integer' []              []; % for backward compatibility
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'subject'     'string'  []              '';
                               'statmode'    'string'  { 'individual' 'common' 'trials' } STUDY.etc.erspparams.statmode}, 'std_erspstatplot');
if isstr(opt), error(opt); end;

% for backward compatibility
% --------------------------
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if ~isempty(opt.comps), 
    opt.subject = STUDY.datasetinfo( STUDY.cluster(opt.clusters).sets(1,opt.comps)).subject;
end;

if ~isempty(opt.subject), statgroup = 'off'; disp('No group statistics for single subject');
else                      statgroup = STUDY.etc.erspparams.statgroup;
end;
if ~isempty(opt.subject), statcond = 'off'; disp('No condition statistics for single subject');
else                      statcond = STUDY.etc.erspparams.statcond;
end;
plotcurveopt = { ...
   'ersplim',    eval( [ 'STUDY.etc.erspparams.' opt.datatype 'lim' ]), ...
   'threshold',  STUDY.etc.erspparams.threshold, ...
   'maskdata',   STUDY.etc.erspparams.maskdata, ...
   'statgroup',  statgroup, ...
   'statcond',   statcond, ...
   'statistics', STUDY.etc.erspparams.statistics };

if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', opt.datatype, 'timerange', opt.timerange, 'statmode', opt.statmode);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', opt.datatype, 'timerange', opt.timerange, 'statmode', opt.statmode);
end;

opt.legend = 'off';
if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);

for index = 1:length(allinds)

    if length(allinds) > 1, subplot(nr,nc,index); end;
    if ~isempty(opt.channels)
        eval( [ 'allersp  = STUDY.changrp(allinds(index)).' opt.datatype 'data;' ]);
        eval( [ 'alltimes = STUDY.changrp(allinds(index)).' opt.datatype 'times;' ]);
        eval( [ 'allfreqs = STUDY.changrp(allinds(index)).' opt.datatype 'freqs;' ]);
    else
        eval( [ 'allersp  = STUDY.cluster(allinds(index)).' opt.datatype 'data;' ]);
        eval( [ 'alltimes = STUDY.cluster(allinds(index)).' opt.datatype 'times;' ]);
        eval( [ 'allfreqs = STUDY.cluster(allinds(index)).' opt.datatype 'freqs;' ]);
    end;

    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(allersp,1)
            for g = 1:size(allersp,2)
                allersp{c,g} = allersp{c,g}(:,subjind);
            end;
        end;
    end;

    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plot({ allfreqs alltimes }, allersp, 'condname', STUDY.condition, 'subject', opt.subject, 'legend', opt.legend, ...
                                      'datatype', opt.datatype,'plotmode', opt.plotmode, 'groupname', STUDY.group, 'plotx', opt.plottf, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
