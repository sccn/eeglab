% pop_importpres() - append presentation file events into eeglab()
%
% Usage:
%   >> EEGOUT = pop_importpres( EEGIN, filename );
%
% Inputs:
%   EEGIN          - input dataset
%   filename       - file name
% 
% Outputs:
%   EEGOUT         - data structure
%
% Note: 1) if they are pre-existing events in the input dataset,
%       this function will recalculate the latency of the events
%       in the presentation file, so that they match the one
%       of the pre-existing events.
%       2) this function call pop_importevent()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2002
%
% See also: eeglab(), pop_importevent()

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
% Revision 1.2  2002/08/07 17:40:24  arno
% header
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, command] = pop_importpres(EEG, filename); 
command = '';

if nargin < 1 
    help pop_importpres;
    return
end;

if nargin < 2 
	% ask user
	[filename, filepath] = uigetfile('*', 'Choose a file from Presentation -- pop_importpres()'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath filename];
end;

EEG = pop_importevent(EEG, 'append', 'no', 'event', filename, 'timeunit', 1E-4, 'skipline', -3, 'align', 0, 'fields', ...
    { 'pres_trial', 'stimulus', 'type', 'latency', 'ttime', 'uncertainty1', 'duration', 'uncertainty2', 'reqtime', 'reqdur' });

command = sprintf('EEG = pop_importpres(%s, ''%s'');', inputname(1), filename); 

return;
