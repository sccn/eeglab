% pop_snapread() - load an EEG SnapMaster file (pop out window if no arguments).
%
% Usage:
%   >> [dat] = pop_snapread( filename, gain);
%
% Graphic interface:
%   "Relative gain" - [edit box] to compute the relative gain, fisrt look at
%                   the text header of the snapmater file with a text editor. 
%                   Find the recording unit, usually in volts (UNITS field).  
%                   Then, find the voltage range in the "CHANNEL.RANGE" [cmin cmax]
%                   field. Finally, determine the gain of the amplifiers (directly
%                   on the machine, not in the header file).
%                   Knowing that the recording precision is 12 bits. The folowing
%                   formula 
%                                    1/2^12*[cmax-cmin]*1e6/gain 
%                   returns the relative gain. You have to compute it and enter
%                   it in the edit box. Enter 1, for preserving the data file units. 
%                   (note that if the voltage range is not the same for all channels
%                   or if the CONVERSION.POLY field in the file header
%                   is not "0 + 1x" for all channels,  you will have to load the data 
%                   using snapread() and scale manually all channels, then import
%                   the Matlab array into EEGLAB).
%
% Inputs:
%   filename       - SnapMaster file name
%   gain           - relative gain. See graphic interface help.
% 
% Outputs:
%   dat            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL/Salk Institute, 13 March 2002
%
% See also: eeglab(), snapread()

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, command] = pop_snapread(filename, gain); 
command = '';
EEG = [];

if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.SMA', 'Choose a SnapMaster file -- pop_snapread()'); 
	if filename == 0 return; end
	filename = [filepath filename];
     
    promptstr    = { 'Relative gain (see help)' };
    inistr       = { '400' };
    result       = inputdlg2( promptstr, 'Import SnapMaster file -- pop_snapread()', 1,  inistr, 'pop_snapread');
    if length(result) == 0 return; end
    gain   = eval( result{1} );

end

if exist('gain') ~= 1
    gain = 1;
end

% load datas
% ----------
EEG = eeg_emptyset;
[EEG.data,params,events, head] = snapread(filename);  

EEG.data            = EEG.data*gain;
EEG.comments        = [ 'Original file: ' filename ];
EEG.filepath        = '';
EEG.setname 		= 'SnapMaster file';
EEG.nbchan          = params(1);
EEG.pnts            = params(2);
EEG.trials          = 1;
EEG.srate           = params(3);
EEG.xmin            = 0; 

A = find(events ~= 0);
if ~isempty(A)
    EEG.event = struct( 'type', mattocell(events(A), [1], ones(1,length(events(A)))), ...
                        'latency', mattocell(A(:)', [1], ones(1,length(A))) );
end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');
command = sprintf('EEG = pop_snapread(''%s'', %f);', filename, gain); 

return;
