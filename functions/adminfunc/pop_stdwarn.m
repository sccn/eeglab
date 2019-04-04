% pop_stdwarn() - check memory options and issue warning for studies.
%
% Usage: >> pop_stdwarn;
%
% Author: Arnaud Delorme, CERCO, 2007
%
% See also: eeg_options()

% Copyright (C) Arnaud Delorme, CERCO, 2007, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
    if strcmpi(res, 'Cancel'),return; end
end
