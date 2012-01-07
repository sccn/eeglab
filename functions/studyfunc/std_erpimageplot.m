% std_erpimageplot() - Commandline function to plot cluster ERPimage or channel erpimage.
%
% Usage:    
%   >> [STUDY] = std_erpimageplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY data times freqs pgroup pcond pinter] = ...
%                std_erpimageplot(STUDY, ALLEEG ...);
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets in the STUDY. 
%                Note: ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Additional help:
% Inputs and output of this function are strictly identical to the std_erspplot(). 
% See the help message of this function for more information. 
%
% See also: std_erspplot()
%
% Authors: Arnaud Delorme, UCSD/CERCO, August, 2011-

% Copyright (C) Arnaud Delorme, arno@ucsd.edu
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

function [STUDY, allitc, alltimes, allfreqs, pgroup, pcond, pinter] = std_erpimageplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_erpimageplot;
    return;
end;

[STUDY allitc alltimes allfreqs pgroup pcond pinter ] = std_erspplot(STUDY, ALLEEG, 'datatype', 'erpim', varargin{:});
