% pop_crossf() - Return estimates and plots of event-related spectral coherence 
%
% Usage:
%       >> pop_crossf(EEG, typeproc, num1, num2, tlimits,cycles,
%                                 'key1',value1,'key2',value2, ... );   
% Inputs:            
%   INEEG    - Input EEG dataset
%   typeproc - Type of processing: 
%                  1 = process two raw-data channels,
%                  0 = process two ICA components
%   num1     - First component or channel number
%   num2     - Second component or channel number
%   tlimits  - [mintime maxtime] Sub-epoch time limits in ms
%   cycles   -   >0 -> Number of cycles in each analysis wavelet 
%                 0 -> Use FFTs (with constant window length)
%
% Optional inputs: As for crossf().  See >> help crossf
%    
% Outputs: Same as crossf(). No outputs are returned when a
%          window pops-up to ask for additional arguments
%
% Author: Arnaud Delorme, CNL / Salk Institute, 11 March 2002
%
% See also: timef(), eeglab() 

% Copyright (C) 11 March 2002 arno@salk.edu, Arnaud Delorme, CNL / Salk Institute
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
	popup = ischar(num1) | isempty(num1);
	if ischar(num1)
		lastcom = num1;
	end
end

% pop up window
% -------------
if popup
	[txt vars] = gethelpvar('timef.m');
	
	geometry = { [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [1 0.5 0.5] [0.92 0.1 0.78] [1 0.5 0.5] [1 0.8 0.2] [1] [1 1]};
    uilist = { { 'Style', 'text', 'string', fastif(typeproc, 'First channel number', 'First component number'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,3,[],'1') } {} ...
			   { 'Style', 'text', 'string', fastif(typeproc, 'Second channel number', 'Second component number'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,4,[],'2') } {} ...
			   { 'Style', 'text', 'string', 'Epoch time range [min max] (msec)', 'fontweight', 'bold', ...
				 'tooltipstring', 'Sub epoch time limits' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,5,[],[ int2str(EEG.xmin*1000) ' ' int2str(EEG.xmax*1000) ]) } {} ...
			   { 'Style', 'text', 'string', 'Wavelet cycles (0->FFT, see >> help timef)', 'fontweight', 'bold', ...
				 'tooltipstring', context('cycles',vars,txt) } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,6,[],'3 0.5') } {} ...
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
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotphase','present',0), 'string', ...
				 'Plot coherence phase', 'tooltipstring', ...
				 'Plot coherence phase image in the lower panel' } ...
			 };

	result = inputgui( geometry, uilist, 'pophelp(''pop_crossf'');', ...
					   fastif(typeproc, 'Plot channel cross-coherence -- pop_crossf()', ...
							  'Plot component cross-coherence -- pop_crossf()'));
	if length( result ) == 0 return; end

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
	if isfield(EEG.chanlocs, 'theta')
        if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end
		if typeproc == 1
			options = [options ', ''topovec'', [' int2str([num1 num2]) ...
                       '], ''elocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo' ];
		else % typeproc == 0
			options = [options ', ''topovec'', EEG.icawinv(:, [' int2str([num1 num2]) ...
                       '])'', ''elocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo' ];
		end
	end
    
    % add title
    % ---------
	if isempty( findstr(  'title', result{7}))
        if ~isempty(EEG.chanlocs) && typeproc
            chanlabel1 = EEG.chanlocs(num1).labels;
            chanlabel2 = EEG.chanlocs(num2).labels;
        else
            chanlabel1 = int2str(num1);
            chanlabel2 = int2str(num2);
        end
		if result{5}
            options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') chanlabel1 '-' chanlabel2 ...
					' Coherence'''];
        else
            options = [options ', ''title'',' fastif(typeproc, '''Channel ', '''Component ') chanlabel1 '-' chanlabel2 ...
					' Phase Coherence''' ];
		end
	end
	if ~isempty( result{6} )
		options      = [ options ', ''alpha'',' result{6} ];
	end
	if ~isempty( result{7} )
		  options = [ options ',' result{7} ];
	end
	if ~result{8}
		options = [ options ', ''plotersp'', ''off''' ];
	end
	if ~result{9}
		options = [ options ', ''plotphase'', ''off''' ];
	end
    figure; try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end; 
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
	tmpsig1 = EEG.data(num1,pointrange,:);
	tmpsig2 = EEG.data(num2,pointrange,:);
else
	if ~isempty( EEG.icasphere )
        tmpsig1 = eeg_getdatact(EEG, 'component', num1, 'samples', pointrange);
        tmpsig2 = eeg_getdatact(EEG, 'component', num2, 'samples', pointrange);
	else
		error('You must run ICA first');
	end;	
end;	 
tmpsig1 = reshape( tmpsig1, 1, size(tmpsig1,2)*size(tmpsig1,3));
tmpsig2 = reshape( tmpsig2, 1, size(tmpsig2,2)*size(tmpsig2,3));

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
varargout{1} = sprintf('figure; pop_crossf( EEG, %d, %d, %d, [%s], [%s] %s);', ...
          typeproc, num1, num2, int2str(tlimits), num2str(cycles), options);

%options = [ options ', ''ydir'', ''norm''' ];
com = sprintf( '%s crossf( tmpsig1, tmpsig2, length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles %s);', outstr, options);
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
