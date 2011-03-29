% pop_importev2() - merge a neuroscan EV2 file with input dataset
%                  (pop out window if no arguments).
%
% Usage:
%   >> OUTEEG = pop_importev2( INEEG ); % pop-up window mode
%   >> OUTEEG = pop_importev2( INEEG, filename);
%
% Inputs:
%   INEEG          - input EEGLAB data structure
%   filename       - file name
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2007

% Copyright (C) 2007 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, command] = pop_importev2(EEG, filename); 
command = '';

if nargin < 1
	help  pop_importev2;
	return;
end;	

if nargin < 2 
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose a EV2 file -- pop_importev2'); 
    drawnow;
	if filename == 0 return; end;
    filename = fullfile(filepath, filename);
end;

% find out if there is a line to skip or not
% ------------------------------------------
fid = fopen(filename, 'r');
tmpl = fgetl(fid);
if isempty(findstr('ype', tmpl)), skipline = 0;
else                              skipline = 1;
end;
fclose(fid);

% load datas
% ----------
tmpevent = EEG.event;
try, oldeventlats = [ tmpevent.latency ]; catch, end;
EEG = pop_importevent(EEG, 'fields', { 'num' 'type' 'response' 'acc' 'RT' 'latency'}, ...
                      'skipline', skipline, 'timeunit', 1E-3, 'align', NaN, 'append', 'no', 'event', filename );

tmpevent = EEG.event;
neweventlats = [ tmpevent.latency ];
if ~exist('oldeventlats'), oldeventlats = neweventlats; end;
len = min(min(length(oldeventlats), length(neweventlats)), 10);
if mean(oldeventlats(1:len) - neweventlats(1:len)) > 1
    error('Wrong alignment of ev2 file with data');
end;

command = sprintf('%s = pop_importev2(%s, %s);', inputname(1), inputname(1), filename); 

return;
