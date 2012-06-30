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

% global variables
% ----------------
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