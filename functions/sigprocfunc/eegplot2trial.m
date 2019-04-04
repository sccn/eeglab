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
   		%end
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
	end
end

% perform the conversion
% ----------------------
tmprejelec   = [];
tmpsig 	    = zeros(1,sweeps);
if ~isempty(TMPREJ)
    nbchan = size(TMPREJ,2)-5;
    TMPREJ = sortrows(TMPREJ,1);
    TMPREJ(find(TMPREJ(:,1) == 1),1) = 0;
	tmpsig = TMPREJ(:,1)''/pnts+1;
	I = find(tmpsig == round(tmpsig));
	tmp = tmpsig(I);
	tmpsig = zeros(1,sweeps);
	tmpsig(tmp) = 1;
    I = find(tmpsig);
	tmprejelec   = zeros( nbchan, sweeps);
	tmprejelec(:,I) = TMPREJ(:,6:nbchan+5)';
end

return;
