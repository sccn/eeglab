% eeg_helphelp() - how to use EEGLAB help.
%
% Menus:  =============================================================
% Each EEGLAB menu calls a Matlab function from the commandline. If this 
% function pops up a graphic interface window, the figure title usually 
% contains the name of the function that the window will call. EEGLAB help 
% files are thus simply the collection of all the help files of the called 
% functions. 
%
% Convention for function calls: ======================================
% When the menu calls a function, it uses the EEG dataset as  an
% argument, sometimes with additional parameters. The function then
% pops-up an interactive window to ask for additional parameter values.
% The advantage of this process is that the same function can be called in 
% two ways, either in interactive (pop_) mode or directly from the commandline. 
% This trick allows EEGLAB to build a history of the commands run under an 
% EEGLAB session. (See the h() function for details). This allows users to 
% build their own EEGLAB macros by copying and pasting commands from the 
% EEGLAB history (from h()) into their own Matlab script files.
%
% Using the EEGLAB help: =========== ???????? =========================
% When using the menus, it is possible to call up the help message of any 
% function that is called from the menu by opening the /Help/EEGLAB menu window.
% The help message of each function that is called is then displayed.
% You must then make the correspondence between the parameters entered through
% the graphic (pop_) interface and the function arguments. Though this
% may seem complicated at the beginning, with experience this  should seem 
% straightforward and will also help the user to learn how to use each function 
% directly from the command line or in user-composed Matlab scripts.
% Many EEGLAB functions do not actually process data, in particular the
% 'pop_' functions. You may then look at the help of the function
% which process the data. For instance, the menu item "/Plot/Channel ERP image"
% calls the function pop_erpimage(). This in turn serves as an interface to 
% the function erpimage().

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

help eeg_helphelp
