% pop_crossf() - Return estimates and plots of event-related (log) spectral
%                coherence amplitudes. 
% Usage:
%       >> pop_crossf(EEG, typeproc, num1, num2, tlimits,cycles,
%                        'key1',value1,'key2',value2, ... );   
%     
% Inputs:            
%   INEEG    - Input EEG dataset
%   typeproc - Type of processing: 
%                1 = process two raw-data channels,
%                0 = process ICA components
%   num1     - First component or channel number
%   num2     - Second component or channel number
%   tlimits  - [mintime maxtime] (in ms) Sub-epoch time limits
%   cycles   - >0 -> Number of cycles in each analysis wavelet 
%               0 -> Use FFTs (with constant window length)
%
% Optional inputs: As for crossf().  See >> help crossf
%    
% Outputs: Same as crossf(). No outputs are returned when a
%          window pops-up to ask for additional arguments
%
% Author: Arnaud Delorme, CNL / Salk Institute, 11 March 2002
%
% See also: timef(), eeglab() 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 11 March 2002 arno@salk.edu, Arnaud Delorme, CNL / Salk Institute
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
% Revision 1.12  2002/05/03 00:35:47  arno
% further check
%
% Revision 1.11  2002/04/25 02:54:18  arno
% debugging topovec
%
% Revision 1.10  2002/04/24 22:08:58  arno
% adding tabs
%
% Revision 1.9  2002/04/24 21:53:09  scott
% [same] -sm
%
% Revision 1.8  2002/04/24 21:37:25  scott
% editing topovec call -sm
%
% Revision 1.7  2002/04/24 21:15:40  scott
% editing topovec call -sm
%
% Revision 1.6  2002/04/24 21:13:45  scott
% adjust topovec args -sm
%
% Revision 1.5  2002/04/24 21:11:18  scott
% added topoplots to plot -sm
%
% Revision 1.4  2002/04/10 01:27:12  arno
% padratio=4 default
%
% Revision 1.3  2002/04/09 19:28:43  arno
% modifying argument passing
%
% Revision 1.2  2002/04/05 23:59:36  arno
% correcting figure title.
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 03-18-02 added title -ad & sm
% 04-04-02 added outputs -ad & sm

function varargout = pop_crossf(EEG, typeproc, num1, num2, tlimits, cycles, varargin );

varargout{1} = '';
% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_crossf;
	return;
end;	
lastcom = [];
if nargin < 3
	popup = 1;
else
	popup = isstr(num1) | isempty(num1);
	if isstr(num1)
		lastcom = num1;
	end;
end;

% pop up window
% -------------
if popup
	[txt vars] = gethelpvar('timef.m');
	
	geometry = { [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [0.92 0.1 0.78] [1 0.5 0.5] [1 0.8 0.2] [1] [1 1]};
    uilist = { { 'Style', 'text', 'string', fastif(typeproc, 'First channel number', 'First component number'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,3,[],'1') } {} ...
			   { 'Style', 'text', 'string', fastif(typeproc, 'Second channel number', 'First component number'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,4,[],'2') } {} ...
			   { 'Style', 'text', 'string', 'Epoch time range [min max] (msec)', 'fontweight', 'bold', ...
				 'tooltipstring', 'Sub epoch time limits' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,5,[],[ int2str(EEG.xmin*1000) ' ' int2str(EEG.xmax*1000) ]) } {} ...
			   { 'Style', 'text', 'string', 'Wavelet cycles (0->FFT, see >> help timef)', 'fontweight', 'bold', ...
				 'tooltipstring', context('cycles',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,6,[],'3') } {} ...
			   { 'Style', 'text', 'string',  '[set]->Linear coher / [unset]->Phase coher', 'fontweight', 'bold', ...
				 'tooltipstring', ['Compute linear inter-trial coherence (coher)' 10 ...
					'OR Amplitude-normalized inter-trial phase coherence (phasecoher)'] } ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'phasecoher','present',1) } { } ...
			   { 'Style', 'text', 'string', 'Bootstrap significance level (Ex: 0.01 -> 1%)', 'fontweight', 'bold', ...
				 'tooltipstring', context('alpha',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,'alpha') } {} ...
			   { 'Style', 'text', 'string', 'Optional timef() arguments (see Help)', 'fontweight', 'bold', ...
				 'tooltipstring', 'See crossf() help via the Help button on the right...' } ...
			   { 'Style', 'edit', 'string', '''padratio'', 4' } ...
			   { 'Style', 'pushbutton', 'string', 'Help', 'callback',  'pophelp(''crossf'');' } ...
			   {} ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotersp','present',0), 'string', ...
				 'Plot coherence amplitude', 'tooltipstring', ...
				 'Plot coherence ampltitude image in the upper panel' } ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotitc','present',0), 'string', ...
				 'Plot coherence phase', 'tooltipstring', ...
				 'Plot coherence phase image in the lower panel' } ...
			 };

	result = inputgui( geometry, uilist, 'pophelp(''pop_crossf'');', ...
					   fastif(typeproc, 'Plot channel cross-coherence -- pop_crossf()', ...
							  'Plot component cross-coherence -- pop_crossf()'));
	if length( result ) == 0 return; end;

	num1     = eval( [ '[' result{1} ']' ] ); 
	num2     = eval( [ '[' result{2} ']' ] ); 
	tlimits	 = eval( [ '[' result{3} ']' ] ); 
	cycles	 = eval( [ '[' result{4} ']' ] );
    if result{5}
    	options = [',''type'', ''coher''' ];
    else
		options = [',''type'', ''phasecoher''' ];
    end;	
	
    % add topoplot
    % ------------
	if ~isempty(EEG.chanlocs)
		if typeproc == 1
			options = [options ', ''topovec'', [' int2str([num1 num2]) '], ''elocs'', EEG.chanlocs' ];
		else % typeproc == 0
			options = [options ', ''topovec'', EEG.icawinv(:, [' int2str([num1 num2]) '])'', ''elocs'', EEG.chanlocs' ];
		end;
	end;
    
    % add title
    % ---------
	if isempty( strmatch(  '''title''', result{7}))
		switch lower(result{5})
		 case 'coher', options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') int2str(num1) '-' int2str(num2) ...
					' Coherence' fastif(~isempty(EEG.setname), [' (' EEG.setname ')'''], '''')];
		 otherwise, options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') int2str(num1) '-' int2str(num2) ...
					' Phase Coherence' fastif(~isempty(EEG.setname), [' (' EEG.setname ')'''],'''') ];
		end;
	end;
	if ~isempty( result{6} )
		options      = [ options ', ''alpha'',' result{6} ];
	end;
	if ~isempty( result{7} )
		  options = [ options ',' result{7} ];
	end;
	if ~result{8}
		options = [ options ', ''plotersp'', ''off''' ];
	end;
	if ~result{9}
		options = [ options ', ''plotitc'', ''off''' ];
	end;
	figure;
else
	options = [];
	for i=1:length( varargin )
		if isstr( varargin{ i } )
			options = [ options ', ''' varargin{i} '''' ];
		else
			if ~isnumeric(varargin) | (size(varargin{i},1) > 1 & size(varargin{i},2) > 1)
				options = [ options ', varargin{' int2str(i) '}' ];
			else
				options = [ options ', [' num2str(varargin{i}) ']' ];	
			end;
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
	tmpsig1 = EEG.data(num1,pointrange,:);
	tmpsig2 = EEG.data(num2,pointrange,:);
else
	if ~isempty( EEG.icasphere )
        eeg_options; % changed from eeglaboptions 3/30/02 -sm
 	    if option_computeica  
    		tmpsig1 = EEG.icaact(num1,pointrange,:);
    		tmpsig2 = EEG.icaact(num2,pointrange,:);
 	    else
            tmpsig1 = (EEG.icaweights(num1,:)*EEG.icasphere)*reshape(EEG.data(:,pointrange,:), EEG.nbchan, EEG.trials*length(pointrange));
            tmpsig2 = (EEG.icaweights(num2,:)*EEG.icasphere)*reshape(EEG.data(:,pointrange,:), EEG.nbchan, EEG.trials*length(pointrange));
		end;
	else
		error('You must run ICA first');
	end;	
end;	 
tmpsig1 = reshape( tmpsig1, 1, size(tmpsig1,2)*size(tmpsig1,3));
tmpsig2 = reshape( tmpsig2, 1, size(tmpsig2,2)*size(tmpsig2,3));

%
% outputs
% -------
outstr = '';
if ~popup
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the datas and generate output command
% --------------------------------------------
if length( options ) < 2
    options = '';
end;
varargout{1} = sprintf('figure; pop_crossf( %s, %d, %d, %d, [%s], %d %s);', ...
          inputname(1), typeproc, num1, num2, int2str(tlimits), cycles, options);
com = sprintf( '%s crossf( tmpsig1, tmpsig2, length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles %s);',... 
          outstr, options);
eval(com)

return;

% get contextual help
% -------------------
function txt = context(var, allvars, alltext);
	loc = strmatch( var, allvars);
	if ~isempty(loc)
		txt= alltext{loc(1)};
	else
		disp([ 'warning: variable ''' var ''' not found']);
		txt = '';
	end;
