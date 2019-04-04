% pop_averef() - Convert an EEG dataset to average reference.
%                This function is obsolete. See pop_reref() instead.
%
% Usage:
%       >> EEGOUT = pop_averef( EEG );
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 March 2002
%
% See also: eeglab(), reref(), averef()

% Copyright (C) 22 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = pop_averef( EEG, confirm);

[EEG, com] = pop_reref(EEG, []);
return;

com = '';
if nargin < 1
   help pop_averef;
   return;
end;   
if isempty(EEG.data)
    error('Pop_averef: cannot process empty data');
end

if nargin < 2 || confirm == 1
    % which set to save
	% -----------------
	 ButtonName=questdlg2( strvcat('Convert the data to average reference?', ...
								   'Note: ICA activations will also be converted if they exist...'), ...
	        'Average reference confirmation -- pop_averef()', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', return;
	 end
	 confirm = 0;
end

EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
if ~isempty(EEG.icaweights)
	disp('pop_averef(): converting ICA weight matrix to average reference (see >> help averef)');
	[EEG.data EEG.icaweights EEG.icasphere EEG.rmave] = averef(EEG.data,EEG.icaweights,EEG.icasphere);
	EEG.icawinv = [];
	if size(EEG.icaweights,1) > EEG.nbchan
		disp('Warning: one or more channels may have been removed; component weight re-referencing may be inaccurate'); 
	end
	if size(EEG.icasphere,1) <  EEG.nbchan
		disp('Warning: one or more components may have been removed; component weight re-referencing could be inaccurate'); 
	end
else
	EEG.data = averef(EEG.data);
end;	
EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
EEG.averef = 'Yes';
EEG.icaact = [];
EEG = eeg_checkset(EEG);

com = sprintf('EEG = pop_averef( EEG, %d);', confirm);
return;
