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
% Revision 1.4  2002/08/09 01:43:18  arno
% [6~[6~same
%
% Revision 1.3  2002/08/09 01:42:40  arno
% debugging filter over smal time period
%
% Revision 1.2  2002/08/09 00:41:22  arno
% updating for boundaries
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_eegfilt( EEG, locutoff, hicutoff, filtorder);

com = '';
if nargin < 1
	help pop_eegfilt;
	return;
end;	
if isempty(EEG.data)
    disp('Pop_eegfilt() error: cannot filter empty dataset'); return;
end;    
if nargin < 3
	% which set to save
	% -----------------
   	promptstr = { 'Lower edge of the frequency pass band (Hz) (0 -> lowpass)', ...
   				  'Higher edge of the frequency pass band (Hz) (0 -> highpass)', ...
   				  'Filter length in points (default: see >> help pop_eegfilt)' };
	inistr       = { '0', '0', '' };
   	result       = inputdlg2( promptstr, 'Filter the data -- pop_eegfilt()', 1,  inistr, 'pop_eegfilt');
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
 
options = { EEG.srate, locutoff, hicutoff, EEG.pnts };
if ~isempty( filtorder )
	options = { options{:} filtorder };
end;

if EEG.trials == 1 
	if ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
		boundaries = strmatch('boundary', {EEG.event.type});
		if isempty(boundaries)
			EEG.data = eegfilt( EEG.data, options{:}); 
		else
			options{4} = 0;
			disp('Pop_eegfilt:finding continuous data boundaries');
			boundaries = [0 round(boundaries-0.5) EEG.pnts];
			for n=1:length(boundaries)-1
				try
					EEGdata(:,boundaries(n)+1:boundaries(n+1)) = ...
						eegfilt(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{:});
				catch
					fprintf('\nPop_eegfilt: data portion from point %d to %d is too small, filter can not be applied\n', boundaries(n),boundaries(n+1));
					disp('Pop_eegfilt: Filter being applied over the whole time range (ignoring discontinuities)');
				end;
			end
		end
	else
		EEG.data = eegfilt( EEG.data, options{:});
	end;
else
	EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
	EEG.data = eegfilt( EEG.data, options{:});
	% note: reshape does not reserv new memory while EEG.data(:,:) does
end;	

com = sprintf( '%s = pop_eegfilt( %s, %s, %s, [%s]);', inputname(1), inputname(1), ...
			num2str( locutoff), num2str( hicutoff), num2str( filtorder ) );
return;			

	 
