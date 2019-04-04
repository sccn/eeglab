% eeg_global() - declare global EEGLAB variables. These variables are
%                used only as global by the main function eeglab(),
%                the function pop_rejmenu() and the history eegh() function.
%           
% Global variables:
%  EEG        - structure containing the current dataset
%  ALLEEG     - array of structures containing all the loaded datasets
%  CURRENTSET - index of the current dataset in the ALLEEG array
%  LASTCOM    - most recent command run by EEGLAB
%  ALLCOM     - cell array containing the EEGLAB session command history 
%               See the history (h) function for more details.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% global variables
% ----------------
globalvars = who('global');
if ~isempty(strmatch('ALLCOM', globalvars, 'exact')) || exist('ALLCOM') ~= 1
    global EEG;		% current dataset 
    global ALLEEG;		% all datasets
    global CURRENTSET;	% current set index

    %global W_MAIN;		% main window
    %global H_MAIN;		% main window
    %global EEGMENU;		% main menu
    global ALLCOM;		% all commands (history)
    global LASTCOM;		% last command
    global STUDY;
    global CURRENTSTUDY;
    global PLUGINLIST;

    if exist('eegplugin_erplab.m')
        global ALLERP; % Javier Lopez-Calderon for ERPLAB
    end
end
