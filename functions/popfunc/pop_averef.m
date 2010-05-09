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
end;

if nargin < 2 | confirm == 1
    % which set to save
	% -----------------
	 ButtonName=questdlg2( strvcat('Convert the data to average reference?', ...
								   'Note: ICA activations will also be converted if they exist...'), ...
	        'Average reference confirmation -- pop_averef()', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', return;
	 end;
	 confirm = 0;
end;

EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
if ~isempty(EEG.icaweights)
	disp('pop_averef(): converting ICA weight matrix to average reference (see >> help averef)');
	[EEG.data EEG.icaweights EEG.icasphere EEG.rmave] = averef(EEG.data,EEG.icaweights,EEG.icasphere);
	EEG.icawinv = [];
	if size(EEG.icaweights,1) > EEG.nbchan
		disp('Warning: one or more channels may have been removed; component weight re-referencing may be inaccurate'); 
	end;
	if size(EEG.icasphere,1) <  EEG.nbchan
		disp('Warning: one or more components may have been removed; component weight re-referencing could be inaccurate'); 
	end;
else
	EEG.data = averef(EEG.data);
end;	
EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
EEG.averef = 'Yes';
EEG.icaact = [];
EEG = eeg_checkset(EEG);

com = sprintf('%s = pop_averef( %s, %d);', inputname(1), inputname(1), confirm);
return;
