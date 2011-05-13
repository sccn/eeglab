% eegplot2trial() - convert EEGPLOT rejections into trial and electrode
%                   rejections compatible with EEGLAB format.
%
% Usage:
%   >> [trialrej elecrej] = eegplot2trial( eegplotrej, frames, ...
%                       sweeps, colorin, colorout );
%
% Inputs:
%   eegplotrej - EEGPLOT output (TMPREJ; see eegplot for more details)
%   frames     - number of points per epoch
%   sweeps     - number of trials
%   colorin    - only extract rejection of specific colors (here a n x 3
%                array must be given). Default: extract all rejections.
%   colorout   - do not extract rejection of specified colors.
%
% Outputs:
%   trialrej   - array of 0 and 1 depicting rejected trials (size sweeps)
%   elecrej    - array of 0 and 1 depicting rejected electrodes in 
%                all trials (size nbelectrodes x sweeps )
%
% Note: if colorin is used, colorout is ignored
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegplot(), eeg_multieegplot(), eegplot2event(), eeglab()

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

function [tmpsig, tmprejelec] = eegplot2trial( TMPREJINIT, pnts, sweeps, color, colorout );

if nargin < 3
   help eegplot2trial;
   return;
end;   

% only take into account specific colors
% --------------------------------------
TMPREJINIT(:,3:5) = round(TMPREJINIT(:,3:5)*100)/100;
TMPREJ = [];
if exist('color') == 1 && ~isempty(color) 
	color    = round(color*100)/100;
	for index = 1:size(color,1)
   		tmpcol1 = TMPREJINIT(:,3) + 255*TMPREJINIT(:,4) + 255*255*TMPREJINIT(:,5);
   		tmpcol2 = color(index,1)+255*color(index,2)+255*255*color(index,3);   
   		I = find( tmpcol1 == tmpcol2);   
   		%if isempty(I)
      	%	fprintf('Warning: color [%d %d %d] not found\n', ...
      	%		color(index,1), color(index,2), color(index,3));
   		%end;
   		TMPREJ = [ TMPREJ; TMPREJINIT(I,:)];
   	end;	
else 
	TMPREJ = TMPREJINIT;
	% remove other specific colors
	% ----------------------------
	if exist('colorout') == 1
		colorout = round(colorout*100)/100;
		for index = 1:size(colorout,1)
			tmpcol1 = TMPREJ(:,3) + 255*TMPREJ(:,4) + 255*255*TMPREJ(:,5);
			tmpcol2 = colorout(index,1)+255*colorout(index,2)+255*255*colorout(index,3);   
			I = find( tmpcol1 == tmpcol2);   
			TMPREJ(I,:) = [];   
		end;	
	end;
end;

% perform the conversion
% ----------------------
tmprejelec   = [];
tmpsig 	    = zeros(1,sweeps);
if ~isempty(TMPREJ)
    nbchan = size(TMPREJ,2)-5;
    TMPREJ(find(TMPREJ(:,1) == 1),1) = 0;
	tmpsig = TMPREJ(:,1)''/pnts+1;
	I = find(tmpsig == round(tmpsig));
	tmp = tmpsig(I);
	tmpsig = zeros(1,sweeps);
	tmpsig(tmp) = 1;
    I = find(tmpsig);
	tmprejelec   = zeros( nbchan, sweeps);
	tmprejelec(:,I) = TMPREJ(:,6:nbchan+5)';
end;

return;
