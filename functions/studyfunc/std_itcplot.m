% std_itcplot() - Commandline function to plot cluster ITCs. Either displays mean cluster 
%                 ITCs, or else all cluster component ITCs, plus the mean cluster ITC, in 
%                 one figure per cluster and condition. ITCs can be visualized only if 
%                 component ITCs were calculated and saved in the STUDY EEG datasets.
%                 These can be computed during pre-clustering using the gui-based function
%                 pop_preclust(), or via the equivalent commandline functions 
%                 eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%              >> [STUDY] = std_itcplot(STUDY, ALLEEG, key1, val1, key2, val2);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets in the STUDY. 
%                Note: ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [numeric vector|'all'] -> indices of clusters to plot.
%                'all' -> plot all clusters in STUDY  {default: 'all'}.
%   'comps'    - [numeric vector|'all']  -> indices of the cluster components to plot.
%                'all' -> plot all components in the cluster {default: 'all'}.
%   'mode'     - ['together'|'apart'] plotting mode. 'together' -> mean ITCs 
%                of the clusters are plotted in the same figure,  one per condition. 
%                'apart' -> component ITCs for each cluster are plotted in a separate 
%                figure (per condition) plus the cluster mean ITC. Note this option is 
%                irrelevant if component indices are provided as input.
%                {default: ' together'}. 
%   'figure'   - ['on'|'off'] 'on' -> plot on a new figure; 'off' -> plot in current
%                figure. 'figure','off' is optional for one cluster in 'centroid' mode.
%                Useful for incorporating cluster ITCs into a complex figure.
%                If multiple conditions, only the first condition is displayed, but
%                clicking on the figure will open a new figure with all conditions 
%                plotted. {default: 'on'}. 
% Outputs:
%   STUDY      - the input STUDY set structure modified with plotted cluster 
%                mean ITCs to allow quick replotting (unless cluster means 
%                already exists in the STUDY).  
%   Example:
%           >> [STUDY] = std_itcplot(STUDY,ALLEEG, 'clusters', 'all', 'mode', 'together');
%              % Plot the mean ITCs of all the clusters in STUDY on the same figure. 
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
% Revision 1.25  2006/09/12 18:51:28  arno
% reprogram from scratch (statistics...), backward compatible
%
                            
function [STUDY allitc alltimes ] = std_itcplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_itcplot;
    return;
end;

[STUDY allitc alltimes ] = std_erspplot(STUDY, ALLEEG, 'datatype', 'itc', varargin{:});
