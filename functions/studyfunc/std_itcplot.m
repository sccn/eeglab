% std_itcplot() - Commandline function to plot cluster ITCs. Either displays mean cluster 
%                 ITCs, or else all cluster component ITCs, plus the mean cluster ITC, in 
%                 one figure per cluster and condition. ITCs can be visualized only if 
%                 component ITCs were calculated and saved in the STUDY EEG datasets.
%                 These can be computed during pre-clustering using the gui-based function
%                 pop_preclust(), or via the equivalent commandline functions 
%                 eeg_createdata() and eeg_preclust(). Called by pop_clustedit().
% Usage:    
%   >> [STUDY] = std_itcplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY itcdata itctimes itcfreqs pgroup pcond pinter] = ...
%                std_itcplot(STUDY, ALLEEG ...);
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets in the STUDY. 
%                Note: ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Additional help:
% Inputs and output of this function are strictly identical to the std_erspplot(). 
% See the help message of this function for more information. std_itcplot()
% plots the ITC while std_erspplot() plots the ERSP.
%
% See also: std_erspplot(), pop_clustedit(), pop_preclust()
%
% Authors: Arnaud Delorme, CERCO, August, 2006-

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

function [STUDY, allitc, alltimes, allfreqs, pgroup, pcond, pinter] = std_itcplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_itcplot;
    return;
end;

[STUDY allitc alltimes allfreqs pgroup pcond pinter ] = std_erspplot(STUDY, ALLEEG, 'datatype', 'itc', varargin{:});
