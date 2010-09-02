% pop_rejtrend() - Measure linear trends in EEG data; reject data epochs 
%                  containing strong trends.
% Usage:
%   >> pop_rejtrend( INEEG, typerej); % pop up an interactive window
%   >> OUTEEG = pop_rejtrend( INEEG, typerej, elec_comp, ...
%                winsize, maxslope, minR, superpose, reject,calldisp);
%
% Pop-up window interface:
%   "Electrode|Component number(s)" - [edit box] electrode or component number(s) 
%                 to take into consideration during rejection. Sets the 'elec_comp'
%                 parameter in the command line call (see below).
%   "Slope window width" - [edit box] integer number of consecutive data
%                 points to use in detecting linear trends. Sets the 'winsize' 
%                 parameter in the command line call.
%   "Maximum slope to allow" - [edit box] maximal absolute slope of the 
%                 linear trend to allow in the data. If electrode data, uV/epoch;
%                 if component data, std. dev./epoch. Sets the 'maxslope'
%                 parameter in the command line call.
%   "R-square limit" -[edit box] maximal regression R-square (0 to 1) value
%                 to allow.  Sets the 'minR' parameter in the command line call.
%                 This represents how "line-like" the rejected data should be; 0
%                 accepts everything that meets the slope requirement, 0.9 is visibly
%                 flat.
%   "Display previous rejection marks?" - [edit box] either YES or NO. 
%                 Sets the command line input option 'superpose'.
%   "Reject marked trials?" - [edit box] either YES or NO. 
%                 Sets the command line input option 'reject'.
% Command line inputs:
%   INEEG      - input EEG dataset
%   typerej    - [1|0] data to reject on: 0 = component activations; 
%                1 = electrode data. {Default: 1}.
%   elec_comp  - [e1 e2 ...] electrode|component number(s) to take into 
%                consideration during rejection
%   winsize    - (integer) number of consecutive points
%                to use in detecing linear trends
%   maxslope   - maximal absolute slope of the linear trend to allow in the data
%   minR       - minimal linear regression R-square value to allow in the data
%                (= coefficient of determination, between 0 and 1)
%   superpose  - [0|1] 0 = Do not superpose marks on previous marks
%                stored in the dataset; 1 = Show both types of marks using 
%                different colors. {Default: 0}
%   reject     - [1|0] 0 = Do not reject marked trials but store the 
%                labels; 1 = Reject marked trials. {Default: 1}
%   calldisp   - [0|1] 1 = Open scroll window indicating rejected trials
%                0 = Do not open scroll window.  {Default: 1}
%
% Outputs:
%   OUTEEG     - output dataset with rejected trials marked for rejection
%     Note: When eegplot() is called, modifications are applied to the current 
%     dataset at the end of the call to eegplot() (e.g., when the user presses 
%     the 'Reject' button).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: rejtrend(), eeglab(), eegplot(), pop_rejepoch() 

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

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-07-02 add the eeglab options -ad

function [EEG, com] = pop_rejtrend( EEG, icacomp, elecrange, winsize, ...
				minslope, minstd, superpose, reject, calldisp);

com = '';
if nargin < 1
   help pop_rejtrend;
   return;
end;  
if nargin < 2
   icacomp = 1;
end;  
if icacomp == 0
	if isempty( EEG.icasphere )
	    ButtonName=questdlg( 'Do you want to run ICA now ?', ...
                         'Confirmation', 'NO', 'YES', 'YES');
    	switch ButtonName,
        	case 'NO', disp('Operation cancelled'); return;   
        	case 'YES', [ EEG com ] = pop_runica(EEG);
    	end % switch
	end;
end;	
if exist('reject') ~= 1
    reject = 1;
end;
if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component number(s), Ex: 2 4 5):', ...
                                           'Electrode number(s), Ex: 2 4 5):'), ...
	                'Slope window width (in points)', ... 
			     fastif(icacomp==0, ...
                                'Maximum slope to allow (uV/epoch)', ...
                                'Maximum slope to allow (std. dev./epoch)'), ...
				'R-square limit to allow (0 to 1, Ex: 0.8):', ...
               		'Display previous rejection marks? (YES or NO)', ...
         			'Reject marked trial(s)? (YES or NO)' };
	inistr      = { ['1:' int2str(EEG.nbchan)], ...
					int2str(EEG.pnts),  ...
					'0.5', ...
					'0.3', ...
               		'NO', ...
            		'NO' };

	result       = inputdlg2( promptstr, fastif(~icacomp, 'Trend rejection in component(s) -- pop_rejtrend()', ...
											   'Data trend rejection -- pop_rejtrend()'), 1,  inistr, 'pop_rejtrend');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	winsize      = result{2};
	minslope     = result{3};
	minstd       = result{4};
    calldisp     = 1;
	switch lower(result{5}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{6}), case 'yes', reject=1; otherwise, reject=0; end;
end;
if nargin < 7
    superpose = 0;
    reject = 0;
    calldisp = 1;
end;
if nargin < 9
    calldisp = 1;
end

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	winsize   = eval( [ '[' winsize ']' ]  );
	minslope  = eval( [ '[' minslope ']' ]  );
	minstd    = eval( [ '[' minstd ']' ]  );
end;

fprintf('Selecting trials...\n');
if icacomp == 1
	[rej tmprejE] = rejtrend( EEG.data(elecrange, :, :), winsize, minslope, minstd);
    rejE = zeros(EEG.nbchan, length(rej));
    rejE(elecrange,:) = tmprejE;
else
    % test if ICA was computed or if one has to compute on line
    % ---------------------------------------------------------
    icaacttmp = eeg_getdatact(EEG, 'component', elecrange); 
    [rej tmprejE] = rejtrend( icaacttmp, winsize, minslope, minstd);
    rejE = zeros(size(icaacttmp,1), length(rej));
    rejE(elecrange,:) = tmprejE;
end
rejtrials = find(rej > 0);
fprintf('%d channel(s) selected\n', size(elecrange(:), 1));
fprintf('%d/%d trial(s) marked for rejection\n', length(rejtrials), EEG.trials);
fprintf('The following trials have been marked for rejection\n');
fprintf([num2str(rejtrials) '\n']);

if calldisp
    if icacomp == 1 macrorej  = 'EEG.reject.rejconst';
        			macrorejE = 'EEG.reject.rejconstE';
    else			macrorej  = 'EEG.reject.icarejconst';
        			macrorejE = 'EEG.reject.icarejconstE';
    end;
	colrej = EEG.reject.rejconstcol;
	eeg_rejmacro; % script macro for generating command and old rejection arrays

	if icacomp == 1
		eegplot( EEG.data(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	else
		eegplot( icaacttmp, 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	end;	
end;
if ~isempty(rej)
	if icacomp	== 1
		EEG.reject.rejconst = rej;
		EEG.reject.rejconstE = rejE;
	else
		EEG.reject.icarejconst = rej;
		EEG.reject.icarejconstE = rejE;
	end;
    if reject
        EEG = pop_rejepoch(EEG, rej, 0);
    end;
end;

%com = sprintf('Indexes = pop_rejtrend( %s, %d, [%s], %s, %s, %s, %d, %d);', ...
%   inputname(1), icacomp, num2str(elecrange),  num2str(winsize), num2str(minslope), num2str(minstd), superpose, reject ); 
com = [ com sprintf('%s = pop_rejtrend(%s,%s);', inputname(1), ...
		inputname(1), vararg2str({icacomp,elecrange,winsize,minslope,minstd,superpose,reject})) ]; 

return;
