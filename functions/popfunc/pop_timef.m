% pop_timef() - Returns estimates and plots of event-related (log) spectral
%           perturbation (ERSP) and inter-trial coherence (ITC) changes 
%           timelocked to a set of input events in one data channel. 
%
% Usage:
%   > pop_timef(EEG, typeproc, num, tlimits,cycles,
%                        'key1',value1,'key2',value2, ... );   
%     
% Inputs:            
%   INEEG    - input EEG dataset
%   typeproc - type of processing. 1 process the raw
%              data and 0 the ICA components
%   num      - component or channel number
%   tlimits  - [mintime maxtime] (ms) sub-epoch time limits
%   cycles   -  >0 -> Number of cycles in each analysis wavelet 
%               0 -> Use FFTs (with constant window length)
%
% Optional inputs:
%    See the timef() function.
%    
% Outputs: same as timef(), no outputs are returned when a
%          window pops-up to ask for additional arguments
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: timef(), eeglab() 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 arno@salk.edu, Arnaud Delorme, CNL / Salk Institute
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
% Revision 1.3  2002/04/06 03:43:36  arno
% adding topoplot options
%
% Revision 1.2  2002/04/05 23:59:06  arno
% correcting figure title
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-08-02 add eeglab option & optimize variable sizes -ad
% 03-10-02 change timef call -ad
% 03-18-02 added title -ad & sm
% 04-04-02 added outputs -ad & sm

function varargout = pop_timef(EEG, typeproc, num, tlimits, cycles, varargin );

varargout{1} = '';
% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_timef;
	return;
end;	

% pop up window
% -------------
if nargin < 3
	promptstr    = { fastif(typeproc, 'Channel number:', 'Component number:') ...
					 'Epoch time range [min max] (msec)' ...
	                 'Wavelet cycles (0->FFT, see >> help timef)' ...
	                 ['Compute linear inter-trial coherence (coher)' 10 'OR inter-trial phase coherence (phasecoher)'] ...
	                 'Bootstrap significance level (Ex: 0.01 -> 1%)' ...
	                 'Optional timef parameters (see >> help timef)' };
	inistr       = {  '1' ...
					  [ int2str(EEG.xmin*1000) ' ' int2str(EEG.xmax*1000) ] ...
					  '0' ...
					  'phasecoher' '' '''padratio'', 4'};
	result       = inputdlg( promptstr, fastif(typeproc, 'Plot channel time frequency -- pop_timef()', ...
	                'Plot component time frequency -- pop_timef()'), 1,  inistr );
	if length( result ) == 0 return; end;

	num	     = eval( [ '[' result{1} ']' ] ); 
	tlimits	 = eval( [ '[' result{2} ']' ] ); 
	cycles	 = eval( [ '[' result{3} ']' ] );
    switch lower(result{4})
    	case 'coher',       options = [',''type'', ''coher''' ];
    	case 'phasecoher',   options = [',''type'', ''phasecoher''' ];
        otherwise, error('Invalid type of coherence');
    end;	
	
    % add topoplot
    % ------------
    if ~isempty(EEG.chanlocs)
      if typeproc == 1
	options = [options ', ''topovec'', int2str(num), ''elocs'', EEG.chanlocs' ];
      else
	options = [options ', ''topovec'', EEG.icawinv(:,int2str(num)), ''elocs'', EEG.chanlocs' ];
      end;
    end;
    
    % add title
    % ---------
	if isempty( strmatch(  '''title''', result{6}))
	    switch lower(result{4})
	       case 'coher', options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') int2str(num) ...
	       ' power and inter-trial coherence' fastif(~ isempty(EEG.setname), [' (' EEG.setname ')''' ], '''') ];
	       otherwise, options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') int2str(num) ...
	       ' power and inter-trial phase coherence' fastif(~ isempty(EEG.setname), [' (' EEG.setname ')''' ], '''') ];
	    end;
	end;
	if ~isempty( result{5} )
		options      = [ options ', ''alpha'',' result{5} ];
	end;
	if ~isempty( result{6} )
		  options = [ options ',' result{6} ];
	else 	
	end;
	figure;
else
	options = [];
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else
			options = [ options ', [' num2str(varargin{i}) ']' ];
		end;
	end;	
end;

% compute epoch limits
% --------------------
if isempty(tlimits)
	tlimits = [EEG.xmin, EEG.xmax];
end;	
pointrange1 = max((tlimits(1)/1000-EEG.xmin)*EEG.srate, 1);
pointrange2 = min((tlimits(2)/1000-EEG.xmin)*EEG.srate, EEG.pnts);
pointrange = [pointrange1:pointrange2];

% call function sample either on raw data or ICA data
% ---------------------------------------------------
if typeproc == 1
	tmpsig = EEG.data(num,pointrange,:);
else
	if ~isempty( EEG.icasphere )
        eeg_options; % changed from eeglaboptions 3/30/02 -sm
 	    if option_computeica  
    		tmpsig = EEG.icaact(num,pointrange,:);
 	    else
            tmpsig = (EEG.icaweights(num,:)*EEG.icasphere)*reshape(EEG.data(:,pointrange,:), EEG.nbchan, EEG.trials*length(pointrange));
        end;
	else
		error('You must run ICA first');
	end;	
end;	 
tmpsig = reshape( tmpsig, length(num), size(tmpsig,2)*size(tmpsig,3));

% outputs
% -------
outstr = '';
if nargin >= 3
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the datas and generate output command
% --------------------------------------------
options
if length( options ) < 2
    options = '';
end;
varargout{1} = sprintf('figure; pop_timef( %s, %d, %d, [%s], %d %s);', inputname(1), typeproc, num, ...
			int2str(tlimits), cycles, options);
com = sprintf('%s timef( tmpsig(:, :), length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles %s);', outstr, options);
eval(com)

return;
