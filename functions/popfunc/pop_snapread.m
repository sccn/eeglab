% pop_snapread() - load an EEG SnapMaster file (pop out window if no arguments).
%
% Usage:
%   >> [dat] = pop_snapread( filename, gain);
%
% Graphic interface:
%   "Relative gain" - [edit box] to compute the relative gain, vizualise
%                   the text header of the snapmater file with a text editor. 
%                   Find the recording unit, usually in volts (UNITS field)).  
%                   Then, find the voltage range in the "CHANNEL.RANGE" [cmin cmax]
%                   field. Finally, look at the gain of the amplifier (directly
%                   on the machine, not in the header file).
%                   Knowing that the recording precision is 12 bits. The folowing
%                   formula 
%                                    1/2^12*[cmax-cmin]*1e6/gain 
%                   returns the relative gain
%                   you have to enter in the edit box. For a gain of 50 and
%                   a channel range of of [-2.5 2.5], the relative gain is 400. 
%                   (note that if the voltage range is not the same for all channels
%                   or if the CONVERSION.POLY field in the file header
%                   is not "0 + 1x" for all channels,  you will have to load the data 
%                   using snapread() and scale manually all channels).
%
% Inputs:
%   filename       - SnapMaster file name
%   gain           - relative gain. See graphic interface help.
% 
% Outputs:
%   dat            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 13 March 2002
%
% See also: eeglab(), snapread()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, command] = pop_snapread(filename, gain); 
command = '';
EEG = [];

if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.SMA', 'Choose a SnapMaster file -- pop_snapread()'); 
	if filename == 0 return; end;
	filename = [filepath filename];
     
    promptstr    = { 'Relative gain (see help)' };
    inistr       = { '400' };
    result       = inputdlg2( promptstr, 'Import SnapMaster file -- pop_snapread()', 1,  inistr, 'pop_snapread');
    if length(result) == 0 return; end;
    gain   = eval( result{1} );

end;

if exist('gain') ~= 1
    gain = 1;
end;

% load datas
% ----------
EEG = eeg_emptyset;
[EEG.data,params,events, head] = snapread(filename);  

EEG.data            = EEG.data*gain;
EEG.filename        = filename;
EEG.filepath        = '';
EEG.setname 		= 'SnapMaster file';
EEG.nbchan          = params(1);
EEG.pnts            = params(2);
EEG.srate           = params(3);
EEG.xmin            = 0; 

A = find(events ~= 0);
if ~isempty(A)
    EEG.event = struct( 'type', mat2cell(events(A), [1], ones(1,length(events(A)))), ...
                    'latency', mat2cell(A(:)', [1], ones(1,length(A))) );
end;

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_snapread(''%s'', %f);', filename, gain); 

return;
