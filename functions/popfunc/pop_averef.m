% pop_averef() - Convert an EEG dataset to average reference.
%
% Usage:
%       >> EEGOUT = pop_averef( EEG, confirm);
%
% Inputs:
%   EEG         - input dataset
%   confirm     - [0|1] ask for confirmation
%
% Inputs:
%   EEGOUT      - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 March 2002
%
% See also: eeglab(), averef()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.2  2002/04/11 17:58:20  arno
% computing average reference of components
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, com] = pop_averef( EEG, confirm);

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
	 ButtonName=questdlg( ...
['Convert the data to average reference?' 10 'Note: ICA activations will also be converted if they exist...'], ...
	        'Average reference confirmation -- pop_averef()', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', return;
	 end;
	 confirm = 0;
end;

EEG.data = reshape(EEG.data, EEG. nbchan, EEG.pnts*EEG.trials);
if ~isempty(EEG.icaweights)
	disp('Pop_averef: converting ICA weight matrix to average reference (see >> help averef)');
	[EEG.data EEG.icaweights EEG.icasphere] = averef(EEG.data,EEG.icaweights,EEG.icasphere);
else
	EEG.data = averef(EEG.data);
end;	
EEG.data = reshape(EEG.data, EEG. nbchan, EEG.pnts, EEG.trials);
EEG.averef = 'Yes';
EEG.icaact = [];
EEG = eeg_checkset(EEG);

com = sprintf('pop_averef( %s, %d);', inputname(1), confirm);
return;
