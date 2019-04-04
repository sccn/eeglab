% eegplot2event() - convert EEGPLOT rejections into events
%                   compatible with the eeg_eegrej function for rejecting 
%                   continuous portions of datasets.
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
%   events     - array of events compatible with the eeg_eegrej function 
%                for rejecting continuous portions of datasets.
%
% Example:
% NEWEEG = eeg_eegrej(EEG,eegplot2event(TMPREJ, -1));
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegplot(), eeg_multieegplot(), eegplot2trial(), eeglab()

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

function events = eegplot2event( TMPREJ, type, color, colorout )

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
   		end
   		TMPREJ = TMPREJ(I,:);   
   	end;	
end

% remove other specific colors
% ----------------------------
if exist('colorout') == 1
	for index = 1:size(colorout,1)
   		tmpcol1 = TMPREJ(:,3) + 255*TMPREJ(:,4) + 255*255*TMPREJ(:,5);
   		tmpcol2 = colorout(index,1)+255*colorout(index,2)+255*255*colorout(index,3);   
   		I = find( tmpcol1 ~= tmpcol2);   
   		TMPREJ = TMPREJ(I,:);   
   	end;	
end

events = [];
if ~isempty(TMPREJ)
   events = TMPREJ(:,1:5);
   events = [ type*ones(size(events,1), 1) ones(size(events,1), 1) round(events(:,1:2)) events(:,3:5)]; 
end

return;
