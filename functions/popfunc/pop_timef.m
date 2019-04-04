% pop_timef() - Returns estimates and plots of event-related (log) spectral
%           perturbation (ERSP) and inter-trial coherence (ITC) changes 
%           timelocked to a set of input events in one data channel. 
%
% Usage:
%   >> pop_timef(EEG, typeplot); % pop_up window
%   >> pop_timef(EEG, typeplot, lastcom); % pop_up window
%   >> pop_timef(EEG, typeplot, channel); % do not pop-up
%   >> pop_timef(EEG, typeproc, num, tlimits,cycles,
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

% Copyright (C) 2002 arno@salk.edu, Arnaud Delorme, CNL / Salk Institute
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
lastcom = [];
if nargin < 3
	popup = 1;
else
	popup = ischar(num) | isempty(num);
	if ischar(num)
		lastcom = num;
	end
end

% pop up window
% -------------
if popup
	[txt vars] = gethelpvar('timef.m');
	
	geometry = { [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [0.92 0.1 0.78] [1 0.5 0.5] [1 0.8 0.2] [1] [1 1]};
    uilist = { ...
                           { 'Style', 'text', 'string', fastif(typeproc, 'Channel number', 'Component number'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,3,[],'1') } {} ...
			   { 'Style', 'text', 'string', 'Epoch time range [min max] (msec)', 'fontweight', 'bold', ...
				 'tooltipstring', 'Sub epoch time limits' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,4,[],[ int2str(EEG.xmin*1000) ' ' int2str(EEG.xmax*1000) ]) } {} ...
			   { 'Style', 'text', 'string', 'Wavelet cycles (0->FFT, see >> help timef)', 'fontweight', 'bold', ...
				 'tooltipstring', context('cycles',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,5,[],'3 0.5') } {} ...
			   { 'Style', 'text', 'string',  '[set]->Linear coher / [unset]->Phase coher', 'fontweight', 'bold', ...
				 'tooltipstring', ['Compute linear inter-trial coherence (coher)' 10 ...
					'OR Amplitude-normalized inter-trial phase coherence (phasecoher)'] } ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'phasecoher','present',1) } { } ...
			   { 'Style', 'text', 'string', 'Bootstrap significance level (Ex: 0.01 -> 1%)', 'fontweight', 'bold', ...
				 'tooltipstring', context('alpha',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,'alpha') } {} ...
			   { 'Style', 'text', 'string', 'Optional timef() arguments (see Help)', 'fontweight', 'bold', ...
				 'tooltipstring', 'See timef() help via the Help button on the right...' } ...
			   { 'Style', 'edit', 'string', '''padratio'', 4, ''plotphase'',''off''' } ...
			   { 'Style', 'pushbutton', 'string', 'Help', 'callback',  'pophelp(''timef'');' } ...
			   {} ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotersp','present',0), 'string', ...
				 'Plot Event Related Spectral Power', 'tooltipstring', ...
				 'Plot log spectral perturbation image in the upper panel' } ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotitc','present',0), 'string', ...
				 'Plot Inter Trial Coherence', 'tooltipstring', ...
				 'Plot the inter-trial coherence image in the lower panel' } ...
			 };
      % { 'Style', 'edit', 'string', '''padratio'', 4, ''plotphase'', ''off''' } ...
			   %{ 'Style', 'text', 'string',  '[set] -> Plot ITC phase sign', 'fontweight', 'bold', ...
			%	 'tooltipstring', ['Plot the sign (+/-) of inter-trial coherence phase' 10 ...
			%		'as red (+) or blue (-)'] } ...
			%   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotphase','present',1) } { } ...

	result = inputgui( geometry, uilist, 'pophelp(''pop_timef'');', ...
					   fastif(typeproc, 'Plot channel time frequency -- pop_timef()', ...
							  'Plot component time frequency -- pop_timef()'));
	if length( result ) == 0 return; end

	num	     = eval( [ '[' result{1} ']' ] ); 
	tlimits	 = eval( [ '[' result{2} ']' ] ); 
	cycles	 = eval( [ '[' result{3} ']' ] );
    if result{4}
    	options = [ ',''type'', ''coher''' ];
    else
		options = [',''type'', ''phasecoher''' ];
    end;	
	
    % add topoplot
    % ------------
    if isfield(EEG.chanlocs, 'theta')
        if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end
        if typeproc == 1
            options = [options ', ''topovec'', ' int2str(num) ...
                        ', ''elocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo' ];
        else
            options = [options ', ''topovec'', EEG.icawinv(:,' int2str(num) ...
                       '), ''elocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo' ];
        end
    end
    
    % add title
    % ---------
	if isempty( findstr(  '''title''', result{6}))
        if ~isempty(EEG.chanlocs) && typeproc
            chanlabel = EEG.chanlocs(num).labels;
        else
            chanlabel = int2str(num);
        end
	    switch lower(result{4})
	       case 'coher', options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') chanlabel ...
	       ' power and inter-trial coherence' fastif(~ isempty(EEG.setname), [' (' EEG.setname ')''' ], '''') ];
	       otherwise, options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') chanlabel ...
	       ' power and inter-trial phase coherence' fastif(~ isempty(EEG.setname), [' (' EEG.setname ')''' ], '''') ];
	    end
	end
	if ~isempty( result{5} )
		options      = [ options ', ''alpha'',' result{5} ];
	end
	if ~isempty( result{6} )
		  options = [ options ',' result{6} ];
	end
	if ~result{7}
		options = [ options ', ''plotersp'', ''off''' ];
	end
	if ~result{8}
		options = [ options ', ''plotitc'', ''off''' ];
	end
	figure; try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
else
    options = [ ',' vararg2str(varargin) ];
end

% compute epoch limits
% --------------------
if isempty(tlimits)
	tlimits = [EEG.xmin, EEG.xmax]*1000;
end;	
pointrange1 = round(max((tlimits(1)/1000-EEG.xmin)*EEG.srate, 1));
pointrange2 = round(min((tlimits(2)/1000-EEG.xmin)*EEG.srate, EEG.pnts));
pointrange = [pointrange1:pointrange2];

% call function sample either on raw data or ICA data
% ---------------------------------------------------
if typeproc == 1
	tmpsig = EEG.data(num,pointrange,:);
else
	if ~isempty( EEG.icasphere )
        tmpsig = eeg_getdatact(EEG, 'component', num, 'samples', pointrange); 
	else
		error('You must run ICA first');
	end;	
end;	 

% outputs
% -------
outstr = '';
if ~popup
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end
end

% plot the datas and generate output command
% --------------------------------------------
if length( options ) < 2
    options = '';
end
if nargin < 4
    varargout{1} = sprintf('figure; pop_timef( EEG, %d, %d, %s, %s %s);', typeproc, num, ...
                           vararg2str({tlimits}), vararg2str({cycles}), options);
end
com = sprintf('%s timef( tmpsig(:, :), length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles %s);', outstr, options);
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
	end
