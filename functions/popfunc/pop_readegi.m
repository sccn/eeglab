% pop_readegi() - load a EGI EEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_readegi;             % a window pops up
%   >> EEG = pop_readegi( filename );
%
% Inputs:
%   filename       - EGI file name
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 Nov 2002
%
% See also: eeglab(), readegi(), readegihdr()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 12 Nov 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.1  2002/11/13 02:34:22  arno
% Initial revision
%

function [EEG, command] = pop_readegi(filename); 
    
EEG = [];
command = '';
if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.RAW', 'Choose an EGI RAW file -- pop_readegi()'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath filename];
end;

% load datas
% ----------
EEG = eeg_emptyset;

[Head EEG.data Eventdata] = readegi( filename );
if ~isempty(Eventdata) & length(Eventdata) == size(EEG.data,2)
    EEG.data(end+1,:) = Eventdata;
end;
EEG.filename        = filename;
EEG.filepath        = '';
EEG.setname 		= 'EGI file';
EEG.nbchan          = size(EEG.data,1);
EEG.srate           = Head.samp_rate;
EEG.trials          = Head.segments;
EEG.pnts            = Head.segsamps;
EEG.xmin            = 0; 

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_readegi(''%s'');', filename); 

return;
