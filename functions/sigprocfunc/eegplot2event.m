% eegplot2event() - convert EEGPLOT rejections into events
%                   compatible with EEGLAB format for continuous
%                   datasets.
%
% Usage:
%   >> [events] = eegplot2event( eegplotrej, type, colorin, colorout );
%
% Inputs:
%   eegplotrej - EEGPLOT output (TMPREJ; see eegplot for more details)
%   type       - type of the event. Default -1.
%   colorin    - only extract rejection of specific colors (here a n x 3
%                array must be given). Default: extract all rejections.
%   colorout   - do not extract rejection of specified colors.
%
% Outputs:
%   events     - array of events in the EEGLAB format.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegplot(), eeg_multieegplot(), eegplot2trial(), eeglab()

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

function events = eegplot2trial( TMPREJ, type, color, colorout );

if nargin < 1
   help eegplot2event;
   return;
end;   

if nargin < 2
   type = -1;
end;   

% only take into account specific colors
% --------------------------------------
if exist('color') == 1
	for index = 1:size(color,1)
   		tmpcol1 = TMPREJ(:,3) + 255*TMPREJ(:,4) + 255*255*TMPREJ(:,5);
   		tmpcol2 = color(index,1)+255*color(index,2)+255*255*color(index,3);   
   		I = find( tmpcol1 == tmpcol2);   
   		if isempty(I)
      		fprintf('Warning: color [%d %d %d] not found\n', ...
      			color(index,1), color(index,2), color(index,3));
   		end;
   		TMPREJ = TMPREJ(I,:);   
   	end;	
end;

% remove other specific colors
% ----------------------------
if exist('colorout') == 1
	for index = 1:size(colorout,1)
   		tmpcol1 = TMPREJ(:,3) + 255*TMPREJ(:,4) + 255*255*TMPREJ(:,5);
   		tmpcol2 = colorout(index,1)+255*colorout(index,2)+255*255*colorout(index,3);   
   		I = find( tmpcol1 ~= tmpcol2);   
   		TMPREJ = TMPREJ(I,:);   
   	end;	
end;

events = [];
if ~isempty(TMPREJ)
   events = TMPREJ(:,1:5);
   events = [ type*ones(size(events,1), 1) ones(size(events,1), 1) round(events(:,1:2)) events(:,3:5)]; 
end;

return;
