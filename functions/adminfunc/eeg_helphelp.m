% eeg_helphelp() - How to use EEGLAB help.
%
% EEGLAB MENU: 
% Each EEGLAB menu item calls a Matlab function from the commandline. If this 
% function pops up a graphic interface window, the figure title usually 
% contains the name of the function that the window will call. EEGLAB help 
% files are thus simply the collection of all the help files of the called 
% functions. 
%
% FUNCTION CALL CONVENTION:
% When the EEGLAB menu calls a function, it takes the EEG dataset as 
% an argument, sometimes with additional parameters. The function then
% pops-up an interactive window asking for additional parameter values.
% The advantage of this process is that the same function can be called in 
% two ways, either in interactive (pop_) mode or directly from the commandline. 
% This trick allows EEGLAB to build a history of the commands run under an 
% EEGLAB session. (See the eegh() function for details). The command history allows 
% users to build their own EEGLAB macros by copying and pasting commands from 
% the EEGLAB history (using eegh() and pop_saveh()) into new Matlab script files.
%
% EEGLAB HELP WINDOWS:
% The help message of any function may be called from from the EEGLAB menu 
% by opening the 'Help > EEGLAB' menu window. The help message of each function 
% is then displayed. Note that many EEGLAB functions do not actually process 
% data (in particular, the 'pop_' functions). To understand their use, look 
% at the help message for the (non-pop) function they call which actually 
% processes the data. For example, the menu item % "Plot > Channel ERP image" 
% calls the function 'pop_erpimage()'. This function in turn serves as a 
% graphic user interface for the computing and plotting function 'erpimage()'.

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

help eeg_helphelp
