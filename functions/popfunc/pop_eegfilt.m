% pop_eegfilt() - filter dataset
%
% Usage:
%   >> eegout = pop_eegfilt( eegin, locutoff, hicutoff, filtorder);
%
% Inputs:
%   eegin     - input dataset
%   locutoff  - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz) {0 -> highpass}
%   filtorder - length of the filter in points {default 3*fix(srate/locutoff)}
%
% Outputs:
%   eegout   - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegfilt(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_eegfilt( EEG, locutoff, hicutoff, filtorder);

com = '';
if isempty(EEG.data)
    disp('Pop_eegfilt() error: cannot filter empty dataset'); return;
end;    
if nargin < 1
	help pop_eegfilt;
	return;
end;	
if nargin < 3
	% which set to save
	% -----------------
   	promptstr = { 'Lower edge of the frequency pass band (Hz) (0 -> lowpass)', ...
   				  'Higher edge of the frequency pass band (Hz) (0 -> highpass)', ...
   				  'Filter length in points (default: see >> help pop_eegfilt)' };
	inistr       = { '0', '0', '' };
   	result       = inputdlg( promptstr, 'Filter the data -- pop_eegfilt()', 1,  inistr);
	if size(result, 1) == 0 return; end;
	locutoff   	 = eval( result{1} );
	hicutoff 	 = eval( result{2} );
	if locutoff == 0 & hicutoff == 0 return; end;
	if isempty( result{3} )
		filtorder = [];
	else
		filtorder    = eval( result{3} );
	end;
elseif nargin < 4
	filtorder = [];
end;	

if EEG.trials == 1
	epochframe = 0;
else
	epochframe = EEG.pnts;
end;	

if isempty( filtorder )
	EEG.data = eegfilt( EEG.data(:,:), EEG.srate, locutoff, hicutoff, epochframe);
else
	EEG.data = eegfilt( EEG.data(:,:), EEG.srate, locutoff, hicutoff, epochframe, filtorder);
end;

com = sprintf( '%s = pop_eegfilt( %s, %s, %s, [%s]);', inputname(1), inputname(1), ...
			num2str( locutoff), num2str( hicutoff), num2str( filtorder ) );
return;			

	 
