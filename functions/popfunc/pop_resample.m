% pop_resample() - resample dataset (pop up window).
%
% Usage:
%   >> [OUTEEG] = pop_resample( INEEG ); % pop up interactive window
%   >> [OUTEEG] = pop_resample( INEEG, freq);
%
% Graphical interface:
%   The edit box entitled "New sampling rate" contains the frequency in
%   Hz for resampling the data. Entering a value in this box  is the same 
%   as providing it in the 'freq' input from the command line.
%
% Inputs:
%   INEEG      - input dataset
%   freq       - frequency to resample (Hz)  
%
% Outputs:
%   OUTEEG     - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: resample(), eeglab()

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
% Revision 1.2  2002/08/12 02:30:30  arno
% [6~[6~inputdlg2
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-08-02 remove ica activity resampling (now set to []) -ad
% 03-08-02 debug call to function help -ad
% 04-05-02 recompute event latencies -ad

function [EEG, command] = pop_resample( EEG, freq); 

command = '';
if nargin < 1
    help pop_resample;
    return;
end;     
if isempty(EEG.data)
    disp('Pop_resample error: cannot resample empty dataset'); return;
end;    

if nargin < 2 

	% popup window parameters
	% -----------------------
	promptstr    = {['New sampling rate']};
	inistr       = { int2str(EEG.srate) };
	result       = inputdlg2( promptstr, 'Resample current dataset -- pop_resample()', 1,  inistr, 'pop_resample');
	if length(result) == 0 return; end;
	freq         = eval( result{1} );

end;

% finding the best ratio
[p,q] = rat(freq/EEG.srate, 0.0001); % not used right now 

% set variable
% ------------
EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
oldpnts   = EEG.pnts;
EEG.pnts  = ceil( EEG.pnts / (EEG.srate/freq) );
fprintf('resampling data %3.0f Hz\n', freq);

% resample for multiple channels
% -------------------------
tmpeeglab = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
for index1 = 1:size(EEG.data,1)
   fprintf('%d ', index1);	
   sigtmp = squeeze (EEG.data(index1,:, :));
	tmpeeglab(index1,:, :) = resample( sigtmp, p, q ); 
end;
fprintf('\n');	
EEG.data = tmpeeglab;

% recompute all event latencies
% -----------------------------
if isfield(EEG.event, 'latency')
    fprintf('resampling event latencies...\n');
    for index1 = 1:length(EEG.event)
        EEG.event(index1).latency = EEG.event(index1).latency * EEG.pnts /oldpnts;
    end;
end;

% resample for multiple channels ica
EEG.icaact = [];

% store dataset
EEG.srate   = freq;
fprintf('resampling finished\n');

EEG.setname = [EEG.setname ' resampled'];

command = sprintf('EEG = pop_resample( %s, %d);', inputname(1), freq);
return;

% resample for non multiple
% -------------------------
	Y  = [1 2];
	for index1 = 1:size(EEG.data,1)
		fprintf('%d\n', index1);	
		for index3 = 1:size(EEG.data,3)
			X = [1:EEG.pnts];
			XX = linspace( 1, EEG.pnts, new_EEG.pnts);
			tmpsig = [ squeeze(EEG.data(index1, :, index3))' squeeze(EEG.data(index1, :, index3))'];
   			[Xi,Yi,Zi] = griddata(Y,X', tmpsig , Y, XX', 'invdist');   % interpolate data
			tmpeeglab(index1,:, index3) = Zi(:,1);
			fprintf('.');	
		end;
		fprintf('\n');	
	end;
