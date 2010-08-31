% pop_stdwarn() - check memory options and issue warning for studies.
%
% Usage: >> pop_stdwarn;
%
% Author: Arnaud Delorme, CERCO, 2007
%
% See also: eeg_options()

% Copyright (C) Arnaud Delorme, CERCO, 2007, arno@salk.edu
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

res = 0;

eeglab_options;
if ~option_storedisk
    vartext = strvcat('Your memory options currently allow to store all datasets in memory (RAM)!', ' ', ...
                   'If your study contains a large number of datasets, you should change the memory', ...
                   'settings to allow EEGLAB to only read the dataset header (cancel next action and', ...
                   'use menu item "File > Memory" - first checkbox to allow at most one dataset at a', ...
                   'time in memory). Otherwise your computer might run out of memory.', ' ', ...
                   'NOTE that this is a REQUIRED step to load the tutorial study since it does not', ...
                   'contain the EEG data.', ' ');
    
    res = questdlg2(vartext, 'Study warning', 'Cancel', 'Ok', 'Ok');
    if strcmpi(res, 'Cancel'), break; return; end;
end;
