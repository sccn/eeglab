% pop_newtimef() - Returns estimates and plots of event-related (log) spectral
%           perturbation (ERSP) and inter-trial coherence (ITC) changes 
%           timelocked to a set of input events in one data channel. 
%
% Usage:
%   >> pop_newtimef(EEG, typeplot); % pop_up window
%   >> pop_newtimef(EEG, typeplot, lastcom); % pop_up window
%   >> pop_newtimef(EEG, typeplot, channel); % do not pop-up
%   >> pop_newtimef(EEG, typeproc, num, tlimits,cycles,
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
%    See the newtimef() function.
%    
% Outputs: same as newtimef(), no outputs are returned when a
%          window pops-up to ask for additional arguments
%
% Getting the ERSP and ITC output values:
% Simply look up the history using the eegh function (type eegh).
% Then copy and paste the pop_newtimef command call and manually add output
% (see the newtimef function for a list of outputs). For instance
% [ersp itc powbase times frequencies] = pop_newtimef( EEG, ....);
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: newtimef(), eeglab() 

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

% 01-25-02 reformated help & license -ad 
% 03-08-02 add eeglab option & optimize variable sizes -ad
% 03-10-02 change newtimef call -ad
% 03-18-02 added title -ad & sm
% 04-04-02 added outputs -ad & sm

function varargout = pop_newtimef(EEG, typeproc, num, tlimits, cycles, varargin );

varargout{1} = '';
% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_newtimef;
	return;
end;	
lastcom = [];
if nargin < 3
	popup = 1;
else
	popup = isstr(num) | isempty(num);
	if isstr(num)
		lastcom = num;
	end;
end;

% pop up window
% -------------
if popup
	[txt vars] = gethelpvar('newtimef.m');
	
    g = [1 0.3 0.6 0.34];
	geometry = { g g g g g g g g [1.025 1.27] [1] [1.2 1 1.2]};
    uilist = { ...
               { 'Style', 'text', 'string', fastif(typeproc, 'Channel number', 'Component number'), 'fontweight', 'bold'  } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,3,[],'1') 'tag' 'chan'} {} {} ...
               ...
			   { 'Style', 'text', 'string', 'Sub epoch time limits [min max] (msec)', 'fontweight', 'bold' } ...
			   { 'Style', 'edit', 'string', getkeyval(lastcom,4,[],[ int2str(EEG.xmin*1000) ' ' int2str(EEG.xmax*1000) ]) 'tag' 'tlimits' } ...
               { 'Style', 'popupmenu', 'string', 'Use 50 time points|Use 100 time points|Use 150 time points|Use 200 time points|Use 300 time points|Use 400 time points' 'tag' 'ntimesout' 'value' 4} { } ...
			   ...
			   { 'Style', 'text', 'string', 'Frequency limits [min max] (Hz) or sequence', 'fontweight', 'bold' } ...
			   { 'Style', 'edit', 'string', '' 'tag' 'freqs'  } ...
               { 'Style', 'popupmenu', 'string', 'Use limits, padding 1|Use limits, padding 2|Use limits, padding 4|Use actual freqs.' 'tag' 'nfreqs' }  ...
               { 'Style', 'checkbox', 'string' 'Log spaced' 'value' 0 'tag' 'freqscale' } ...
			   ...
			   { 'Style', 'text', 'string', 'Baseline limits [min max] (msec) (0->pre-stim.)', 'fontweight', 'bold' } ...
			   { 'Style', 'edit', 'string', '0' 'tag' 'baseline' } ...
			   { 'Style', 'popupmenu',  'string', 'Use divisive baseline|Use standard deviation' 'tag' 'basenorm' } ...
               { 'Style', 'checkbox', 'string' 'No baseline' 'tag' 'nobase' } ...
               ...
               { 'Style', 'text', 'string', 'Wavelet cycles [min max/fact] or sequence', 'fontweight', 'bold' } ...
               { 'Style', 'edit', 'string', getkeyval(lastcom,5,[],'3 0.5') 'tag' 'cycle' } ...
               { 'Style', 'popupmenu', 'string', 'Use limits|Use actual seq.' 'tag' 'ncycles' } ...
               { 'Style', 'checkbox', 'string' 'Use FFT' 'value' 0 'tag' 'fft' } ...
			   ...
			   { 'Style', 'text', 'string', 'ERSP color limits [max] (min=-max)', 'fontweight', 'bold' } ...
               { 'Style', 'edit', 'string', '' 'tag' 'erspmax'} ...
               { 'Style', 'checkbox', 'string' 'see log power (set)' 'tag' 'scale' 'value' 1} {} ...
			   ...
			   { 'Style', 'text', 'string', 'ITC color limits [max]', 'fontweight', 'bold' } ...
               { 'Style', 'edit', 'string', '' 'tag' 'itcmax'} ...
               { 'Style', 'checkbox', 'string' 'plot ITC phase (set)' 'tag' 'plotphase' } {} ...
			   ...
			   { 'Style', 'text', 'string', 'Bootstrap significance level (Ex: 0.01 -> 1%)', 'fontweight', 'bold' } ...
               { 'Style', 'edit', 'string', getkeyval(lastcom,'alpha') 'tag' 'alpha'} ...
               { 'Style', 'checkbox', 'string' 'FDR correct (set)' 'tag' 'fdr' } {} ...
			   ...
			   { 'Style', 'text', 'string', 'Optional newtimef() arguments (see Help)', 'fontweight', 'bold', ...
				 'tooltipstring', 'See newtimef() help via the Help button on the right...' } ...
			   { 'Style', 'edit', 'string', '' 'tag' 'options' } ...
			   {} ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotersp','present',0), 'string', ...
				 'Plot Event Related Spectral Power', 'tooltipstring', ...
				 'Plot log spectral perturbation image in the upper panel' 'tag' 'plotersp' } ...
			   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotitc','present',0), 'string', ...
				 'Plot Inter Trial Coherence', 'tooltipstring', ...
				 'Plot the inter-trial coherence image in the lower panel' 'tag' 'plotitc' } ...
			   { 'Style', 'checkbox', 'value', 0, 'string', ...
				 'Plot curve at each frequency' 'tag' 'plotcurve' } ...
			 };
      % { 'Style', 'edit', 'string', '''padratio'', 4, ''plotphase'', ''off''' } ...
			   %{ 'Style', 'text', 'string',  '[set] -> Plot ITC phase sign', 'fontweight', 'bold', ...
			%	 'tooltipstring', ['Plot the sign (+/-) of inter-trial coherence phase' 10 ...
			%		'as red (+) or blue (-)'] } ...
			%   { 'Style', 'checkbox', 'value', ~getkeyval(lastcom,'plotphase','present',1) } { } ...

	[ tmp1 tmp2 strhalt result ] = inputgui( geometry, uilist, 'pophelp(''pop_newtimef'');', ...
					   fastif(typeproc, 'Plot channel time frequency -- pop_newtimef()', ...
							  'Plot component time frequency -- pop_newtimef()'));
	if length( tmp1 ) == 0 return; end;

	if result.fft,      result.cycle = '0'; end;
	if result.nobase,   result.baseline = 'NaN'; end;
    
	num	     = eval( [ '[' result.chan    ']' ] ); 
	tlimits	 = eval( [ '[' result.tlimits ']' ] ); 
	cycles	 = eval( [ '[' result.cycle   ']' ] );
    freqs    = eval( [ '[' result.freqs   ']' ] );
    %result.ncycles == 2 is ignored
    
    % add topoplot
    % ------------
    options = [];
    if isfield(EEG.chanlocs, 'theta')
        if ~isfield(EEG, 'chaninfo'), EEG.chaninfo = []; end;
        if typeproc == 1
            if isempty(EEG.chanlocs), caption = [ 'Channel ' int2str(num) ]; else caption = EEG.chanlocs(num).labels; end;
            options = [options ', ''topovec'', ' int2str(num) ...
                        ', ''elocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo, ''caption'', ''' caption '''' ];
        else
            options = [options ', ''topovec'', EEG.icawinv(:,' int2str(num) ...
                       '), ''elocs'', EEG.chanlocs, ''chaninfo'', EEG.chaninfo, ''caption'', [''IC '', num2str(num)]' ];
      end;
    end;
    
	if ~isempty( result.baseline ),  options = [ options ', ''baseline'',[' result.baseline ']' ]; end;
    if ~isempty( result.alpha ),     options = [ options ', ''alpha'',' result.alpha ];   end;
	if ~isempty( result.options ),   options = [ options ',' result.options ];            end;
	if ~isempty( result.freqs ),     options = [ options ', ''freqs'', [' result.freqs ']'   ]; end;
	if ~isempty( result.erspmax ),   options = [ options ', ''erspmax'', [' result.erspmax ']' ]; end;
	if ~isempty( result.itcmax ),    options = [ options ', ''itcmax'','  result.itcmax ];      end;
	if ~result.plotersp,             options = [ options ', ''plotersp'', ''off''' ];     end;
	if ~result.plotitc,              options = [ options ', ''plotitc'' , ''off''' ];     end;
	if result.plotcurve,             options = [ options ', ''plottype'', ''curve''' ];   end;
	if result.fdr,                   options = [ options ', ''mcorrect'', ''fdr''' ];     end;
	if result.freqscale,             options = [ options ', ''freqscale'', ''log''' ];    end;
	if ~result.plotphase,            options = [ options ', ''plotphase'', ''off''' ];    end;
	if ~result.scale,                options = [ options ', ''scale'', ''abs''' ];        end;
    if result.basenorm == 2,         options = [ options ', ''basenorm'', ''on''' ];      end;
    if result.ntimesout == 1,        options = [ options ', ''ntimesout'', 50' ];         end;
    if result.ntimesout == 2,        options = [ options ', ''ntimesout'', 100' ];        end;
    if result.ntimesout == 3,        options = [ options ', ''ntimesout'', 150' ];        end;
    if result.ntimesout == 5,        options = [ options ', ''ntimesout'', 300' ];        end;
    if result.ntimesout == 6,        options = [ options ', ''ntimesout'', 400' ];        end;
    if result.nfreqs == 1,           options = [ options ', ''padratio'', 1' ];           end;    
    if result.nfreqs == 2,           options = [ options ', ''padratio'', 2' ];           end;    
    if result.nfreqs == 3,           options = [ options ', ''padratio'', 4' ];           end;
    if result.nfreqs == 4,           options = [ options ', ''nfreqs'', ' int2str(length(freqs)) ]; end;
    
    % add title
    % ---------
	if isempty( findstr(  '''title''', result.options))
        if ~isempty(EEG.chanlocs) & typeproc
            chanlabel = EEG.chanlocs(num).labels;
        else
            chanlabel = int2str(num);
        end;
	end;
    
	figure; try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
else
    options = [ ',' vararg2str(varargin) ];
end;

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
        eeglab_options; % changed from eeglaboptions 3/30/02 -sm
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
if ~popup
    for io = 1:nargout, outstr = [outstr 'varargout{' int2str(io) '},' ]; end;
    if ~isempty(outstr), outstr = [ '[' outstr(1:end-1) '] =' ]; end;
end;

% plot the datas and generate output command
% --------------------------------------------
if length( options ) < 2
    options = '';
end;
if nargin < 4
    varargout{1} = sprintf('figure; pop_newtimef( %s, %d, %d, [%s], [%s] %s);', inputname(1), typeproc, num, ...
			int2str(tlimits), num2str(cycles), options);
end;
com = sprintf('%s newtimef( tmpsig(:, :), length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles %s);', outstr, options);
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
