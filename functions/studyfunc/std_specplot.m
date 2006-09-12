% std_specplot() - visualizes component cluster spectra, either mean spectra for 
%                  all requested clusters in the same figure, with spectra for 
%                  different conditions (if any) plotted in different colors, 
%                  or spectra for each specified cluster in a separate figure 
%                  for each condition,  showing the cluster component spectra plus 
%                  the mean cluster spectrum (in bold). The spectra can be 
%                  plotted only if component spectra have been computed and 
%                  saved with the EEG datasets in Matlab files "[datasetname].icaspec" 
%                  using pop_preclust() or std_preclust(). Called by pop_clustedit(). 
%                  Calls std_readspec() and internal function std_plotcompspec()
% Usage:    
%              >> [STUDY] = std_specplot(STUDY, ALLEEG, key1, val1, key2, val2, ...);  
% Inputs:
%   STUDY      - STUDY structure comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - vector of EEG dataset structures for the dataset(s) in the STUDY, 
%                typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [int vector] -> cluster numbers to plot.
%                       'all' -> plot all clusters in STUDY.
%                {default: 'all'}.
%   'comps'    - [int vector] -> indices of cluster components to plot.
%                       'all' -> plot all the components in the cluster 
%                {default: 'all'}.
%   'mode'     - ['centroid'|'comps'] plotting mode. In 'centroid' mode, the average 
%                spectra of the requested clusters are plotted in the same figure, 
%                with spectra for  different conditions (if any) plotted in different 
%                colors. In 'comps' mode, spectra for each specified cluster are 
%                plotted in separate figures (per condition), each containing the
%                cluster component spectra plus the mean cluster spectrum in bold.
%                {default: 'centroid'}. Note that this option is irrelevant when 
%                component indices are provided as input.
%   'figure'   - ['on'|'off'] for the 'centroid' mode option, 'on' plots in a new 
%                figure, while 'off'  plots in the current figure. {default: 'on'}
%
% Outputs:
%   STUDY      - the input STUDY set structure modified with the plotted cluster 
%                mean spectra, to allow quick replotting (unless the cluster means 
%                already exists in the STUDY).  
%
%   Example:
%            >> [STUDY] = std_specplot(STUDY,ALLEEG, 'clusters', 2, 'mode', 'comps');
%               % Plot component spectra for Cluster 2 plus the mean cluster spectrum 
%               % (in bold). 
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
                            
function [STUDY, erspbase] = std_specplot(STUDY, ALLEEG, varargin)

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
if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;
if ~isempty(opt.comps), 
    opt.subject = STUDY.datasetinfo( STUDY.cluster(opt.clusters).sets(1,opt.comps)).subject;
end;

if ~isempty(opt.subject), statgroup = 'off'; disp('No group statistics for single subject');
else                      statgroup = STUDY.etc.specparams.statgroup;
end;
if ~isempty(opt.subject), statcond = 'off'; disp('No condition statistics for single subject');
else                      statcond = STUDY.etc.specparams.statcond;
end;
plotcurveopt = { ...
   'ylim',       STUDY.etc.specparams.ylim, ...
   'threshold',  STUDY.etc.specparams.threshold, ...
   'statgroup',  statgroup, ...
   'statcond',   statcond, ...
   'plotgroup',  STUDY.etc.specparams.plotgroup, ...
   'plotcond',   STUDY.etc.specparams.plotcond, ...
   'statistics', STUDY.etc.specparams.statistics };

if ~isempty(opt.channels)
     [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'spec', 'freqrange', opt.freqrange);
else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'spec', 'freqrange', opt.freqrange);
end;

opt.legend = 'off';
if length(allinds) > 1, figure; opt.plotmode = 'condensed'; end;
nc = ceil(sqrt(length(allinds)));
nr = ceil(length(allinds)/nc);

for index = 1:length(allinds)

    if length(allinds) > 1, subplot(nr,nc,index); end;
    if ~isempty(opt.channels)
        erspbase = STUDY.changrp(allinds(index)).specdata;
        allfreqs = STUDY.changrp(allinds(index)).specfreqs;
    else
        erspbase = STUDY.cluster(allinds(index)).specdata;
        allfreqs = STUDY.cluster(allinds(index)).specfreqs;
    end;

    if ~isempty(opt.subject)
        subjind = strmatch(opt.subject, STUDY.subject);
        for c = 1:size(erspbase,1)
            for g = 1:size(erspbase,2)
                erspbase{c,g} = erspbase{c,g}(:,subjind);
            end;
        end;
    end;

    if index == length(allinds), opt.legend = 'on'; end;
    [pgroup pcond pinter] = std_plot(allfreqs, erspbase, 'condname', STUDY.condition, 'legend', opt.legend, ...
                                      'plotmode', opt.plotmode, 'groupname', STUDY.group, 'plotx', opt.plotfreq, 'unitx', 'Hz', ...
                                      'chanlocs', ALLEEG(1).chanlocs, 'plotsubjects', opt.plotsubjects, plotcurveopt{:});
    if length(allinds) > 1, 
        if isempty(opt.channels), title(sprintf('Cluster %d', allinds(index))); 
        else                      title(sprintf('%s', opt.channels{index}));  
        end;
    end;
end;
