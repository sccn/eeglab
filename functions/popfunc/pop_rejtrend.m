% pop_rejtrend() - Detect linear trends in EEG activity and reject the  
%                  epoched trials based on the accuracy of the linear
%                  fit.
% Usage:
%   >> pop_rejtrend( INEEG, typerej); % pop up interactive window
%   >> OUTEEG = pop_rejtrend( INEEG, typerej, elec_comp, ...
%                winsize, maxslope, minR, superpose, reject);
%
% Graphical interface:
%   "Electrode" - [edit box] electrodes or components (number) to take into 
%                 consideration for rejection. Same as the 'elec_comp'
%                 parameter from the command line.
%   "Consecutive alike values" - [edit box] integer determining the
%                 number of consecutive points for the detection of linear
%                 patterns. Same as the 'winsize' parameter from the
%                 command line.
%   "Maximal slope" - [edit box] maximal absolute slope of the linear 
%                trend of the activity for rejection. Same as the 'maxslope'
%                 parameter from the command line.
%   "R square limit" -[edit box] minimal R^2 (0 to 1). Same as 'minR'
%                 parameter from the command line.
%   "Display with previous rejection" - [edit box] can be either YES or 
%                 NO. This edit box corresponds to the command line input 
%                 option 'superpose'.
%   "Reject marked trials" - [edit box] can be either YES or NO. This edit
%                 box corresponds to the command line input option 'reject'.
%
% Inputs:
%   INEEG      - input dataset
%   typerej    - type of rejection (0 = independent components; 1 = eeg
%                data). Default is 1.
%   elec_comp  - [e1 e2 ...] electrodes or components (number) to take into 
%                consideration for rejection
%   winsize    - integer determining the number of consecutive points
%                for the detection of linear patterns
%   maxslope   - maximal absolute slope of the linear trend of the 
%                activity for rejection
%   minR       - minimal R^2 (coefficient of determination between
%                0 and 1)
%   superpose  - 0=do not superpose pre-labelling with previous
%                pre-labelling (stored in the dataset). 1=consider both
%                pre-labelling (using different colors). Default is 0.
%   reject     - 0=do not reject labelled trials (but still store the 
%                labels. 1=reject labelled trials. Default is 0.
%
% Outputs:
%   OUTEEG     - output dataset with labeled rejected sweeps
%     when eegplot is called, modifications are applied to the current 
%     dataset at the end of the call of eegplot (when the user press the 
%     button 'reject').
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: rejtrend(), eeglab(), eegplot(), pop_rejepoch() 

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
% Revision 1.13  2003/02/18 22:43:06  arno
% same.
%
% Revision 1.12  2003/02/18 22:20:21  arno
% update header for GUI
% /
%
% Revision 1.11  2002/11/13 01:19:42  arno
% debug for command line call
%
% Revision 1.10  2002/11/12 23:56:52  luca
% now saving outputs -ad
%
% Revision 1.9  2002/08/12 21:53:22  arno
% text
%
% Revision 1.8  2002/08/12 02:31:36  arno
% inputdlg2
%
% Revision 1.7  2002/08/07 22:39:38  arno
% same
%
% Revision 1.6  2002/08/07 22:28:38  arno
% editing
%
% Revision 1.5  2002/07/30 23:00:28  arno
% new rejection type
%
% Revision 1.4  2002/07/26 17:32:49  arno
% better output command
%
% Revision 1.3  2002/07/26 17:04:24  arno
% adding message
%
% Revision 1.2  2002/07/26 16:35:57  arno
% chaning icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-07-02 add the eeglab options -ad

function [EEG, com] = pop_rejtrend( EEG, icacomp, elecrange, winsize, ...
				minslope, minstd, superpose, reject, topcommand);

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
if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component (number; ex: 2 4 5):', 'Electrode (number; ex: 2 4 5):'), ...
	                'Consecutive alike values (in data points)', ... 
					'Maximal slope (trend) of the activity (unit/epoch):', ...
					'R square limit (0 to 1, ex: 0.8):', ...
               		'Display with previous rejection', ...
         			'Reject marked trial(s) (YES or NO)' };
	inistr      = { ['1:' int2str(EEG.nbchan)], ...
					int2str(EEG.pnts),  ...
					'0.5', ...
					'0.3', ...
               		'NO', ...
            		'NO' };

	result       = inputdlg2( promptstr, fastif(~icacomp, 'Trend rejection in component -- po_rejtrend()', ...
											   'Trend rejection -- po_rejtrend()'), 1,  inistr, 'pop_rejtrend');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	winsize      = result{2};
	minslope     = result{3};
	minstd       = result{4};
	switch lower(result{5}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{6}), case 'yes', reject=1; otherwise, reject=0; end;
end;
if nargin < 7
    superpose = 0;
    reject = 0;
    topcommand = '';
end;

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	winsize   = eval( [ '[' winsize ']' ]  );
	minslope  = eval( [ '[' minslope ']' ]  );
	minstd    = eval( [ '[' minstd ']' ]  );
else
	calldisp = 0;
end;

fprintf('Selecting trials...\n');
if icacomp == 1
	[rej tmprejE] = rejtrend( EEG.data(elecrange, :, :), winsize, minslope, minstd);
    rejE = zeros(EEG.nbchan, length(rej));
    rejE(elecrange,:) = tmprejE;
else
    % test if ICA was computed or if one has to compute on line
    % ---------------------------------------------------------
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
        icaacttmp = EEG.icaact(elecrange, :, :);
	else
        icaacttmp = EEG.icaweights(elecrange,:)*EEG.icasphere*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        icaacttmp = reshape( icaacttmp, length(elecrange), EEG.pnts, EEG.trials);
    end;
    [rej tmprejE] = rejtrend( icaacttmp, winsize, minslope, minstd);
    rejE = zeros(size(icaacttmp,1), length(rej));
    rejE(elecrange,:) = tmprejE;
end;
fprintf('%d channel selected\n', size(elecrange(:), 1));
fprintf('%d/%d trials marked for rejection\n', length(find(rej > 0)), EEG.trials);

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
