% pop_rejkurt() - rejection of artifact in a dataset using kurtosis 
%                 of activity (i.e. to detect peaky distribution of
%                 activity).
%
% Usage:
%   >> pop_rejkurt( INEEG, typerej) % pop-up interative window mode
%   >> [OUTEEG, locthresh, globthresh, nrej] = ...
%		= pop_rejkurt( INEEG, typerej, elec_comp, ...
%                   locthresh, globthresh, superpose, reject, vistype);
%
% Graphical interface:
%   "Electrode" - [edit box] electrodes or components (number) to take into
%                 consideration for rejection. Same as the 'elec_comp'
%                 parameter from the command line.
%   "Single-channel limit" - [edit box] pkurtosis limit in terms of
%                 standard-dev. Same as 'locthresh' command line
%                 parameter.
%   "All-channel limit" - [edit box] kurtosis limit in terms of
%                 standard-dev (all channel regrouped). Same as 
%                 'globthresh' command line parameter.
%   "Display with previous rejection" - [edit box] can be either YES or
%                 NO. This edit box corresponds to the command line input
%                 option 'superpose'.
%   "Reject marked trials" - [edit box] can be either YES or NO. This edit
%                 box corresponds to the command line input option 'reject'.
%   "visualization type" - [edit box] can be either REJECTRIALS or EEGPLOT.
%                 This edit box corresponds to the command line input
%                 option 'vistype'.
% 
% Inputs:
%   INEEG      - input dataset
%   typerej    - type of rejection (0 = independent components; 1 = eeg
%              data). Default is 1. For independent components, before
%              thresholding, the activity is normalized for each 
%              component.
%   elec_comp  - [e1 e2 ...] electrodes (number) to take into 
%              consideration for rejection
%   locthresh  - activity kurtosis limit in terms of standard-dev.
%   globthresh - global limit (for all channel). Same unit as above.
%   superpose  - 0=do not superpose pre-labelling with previous
%              pre-labelling (stored in the dataset). 1=consider both
%              pre-labelling (using different colors). Default is 0.
%   reject     - 0=do not reject labelled trials (but still store the 
%              labels. 1=reject labelled trials. Default is 0.
%   vistype    - visualization type. 0 is rejepoch() and 1 is eegplot()
%              default is 0.  
%
% Outputs:
%   OUTEEG     - output dataset with updated kurtosis array
%   locthresh  - electrodes probability of activity thresholds in terms
%              of standard-dev.
%   globthresh - global threshold (where all electrode activity are 
%              regrouped).
%   nrej       - number of rejected sweeps
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: rejkurt(), rejstatepoch(), pop_rejepoch(), eegplot(), eeglab()  

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
% Revision 1.13  2002/08/14 00:57:48  arno
% debugging
%
% Revision 1.12  2002/08/14 00:51:47  arno
% debug multi-elec
%
% Revision 1.11  2002/08/12 21:52:30  arno
% same
%
% Revision 1.10  2002/08/12 21:52:02  arno
% text
%
% Revision 1.9  2002/08/12 02:33:29  arno
% inputdlg2
%
% Revision 1.8  2002/08/07 23:19:15  arno
% editing
%
% Revision 1.7  2002/08/07 22:39:14  arno
% same
%
% Revision 1.6  2002/08/07 22:35:56  arno
% editing
%
% Revision 1.5  2002/07/30 22:26:35  arno
% new rejection
%
% Revision 1.4  2002/07/26 18:07:57  arno
% debugging
%
% Revision 1.3  2002/07/26 17:21:05  arno
% return command if 2 outputs
%
% Revision 1.2  2002/07/26 16:38:49  arno
% inverting icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-08-02 add eeglab options -ad

function [EEG, locthresh, globthresh, nrej, com] = pop_rejkurt( EEG, icacomp, elecrange, ...
                       		locthresh, globthresh, superpose, reject, vistype, topcommand);
com = '';
if nargin < 1
   help pop_rejkurt;
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
	promptstr   = { [ fastif(icacomp, 'Electrode', 'Component') ' (number; ex: 2 4 5):' ], ...
					[ fastif(icacomp, 'Single-channel', 'Single-component') ' limit (standard deviation: ex: 3 4 2):'], ...
					[ fastif(icacomp, 'All-channel', 'All-component') ' limit (standard deviation: ex: 3 4 2):'], ...
               		'Display with previous rejection', ...
         			'Actually reject marked trial(s) (YES or NO)', ...
         			'visualization type (REJECTRIALS|EEGPLOT)' };
	inistr      = { fastif(icacomp, ['1:' int2str(EEG.nbchan)], ['1:' int2str(size(EEG.icaweights,1))])...
					fastif(icacomp, '3', '5'),  ...
					fastif(icacomp, '3', '5'), ...
               		'YES', ...
            		'NO', ...
            		'REJECTRIALS' };

	result       = inputdlg2( promptstr, fastif(~icacomp, 'Trial rejection using comp. kurtosis -- pop_rejkurt()', 'Trial rejection using data kurtosis -- pop_rejkurt()'), 1,  inistr, 'pop_rejkurt');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	locthresh    = result{2};
	globthresh   = result{3};
	switch lower(result{4}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{5}), case 'yes', reject=1; otherwise, reject=0; end;
	switch lower(result{6}), case 'rejepoch', vistype=0; otherwise, vistype=1; end;
end;

if ~exist('vistype') vistype = 0; end;
if ~exist('reject') reject = 0; end;
if ~exist('superpose') superpose = 1; end;

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	locthresh = eval( [ '[' locthresh ']' ]  );
	globthresh = eval( [ '[' globthresh ']' ]  );
else
	calldisp = 0;
end;

if isempty(elecrange)
	error('No electrode selectionned');
end;	

% compute the joint probability
% -----------------------------
if icacomp == 1
	fprintf('Computing kurtosis for channels...\n');
    if isempty(EEG.stats.kurtE )
		[ EEG.stats.kurtE rejE ] = rejkurt( EEG.data, locthresh, EEG.stats.kurtE, 1); 
	end;
	[ tmp rejEtmp ] = rejkurt( EEG.data(elecrange, :,:), locthresh, EEG.stats.kurtE(elecrange, :), 1); 
    rejE    = zeros(EEG.nbchan, size(rejEtmp,2));
	rejE(elecrange,:) = rejEtmp;
	
	fprintf('Computing all-channel kurtosis...\n');
	tmpdata = permute(EEG.data, [3 1 2]);
	tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
	[ EEG.stats.kurt rej ] = rejkurt( tmpdata, globthresh, EEG.stats.kurt, 1); 
else
	fprintf('Computing joint probability for components...\n');
    % test if ICA was computed
    % ------------------------
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
 	if option_computeica  
    	icaacttmp = EEG.icaact;
	else
        icaacttmp = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        icaacttmp = reshape( icaacttmp, size(icaacttmp,1), EEG.pnts, EEG.trials);
    end;
    if isempty(EEG.stats.icakurtE )
		[ EEG.stats.icakurtE rejE ] = rejkurt( icaacttmp, locthresh, EEG.stats.icakurtE, 1); 
	end;
	[ tmp rejEtmp ] = rejkurt( icaacttmp(elecrange, :,:), locthresh, EEG.stats.icakurtE(elecrange, :), 1); 
	rejE    = zeros(size(icaacttmp,1), size(rejEtmp,2));
	rejE(elecrange,:) = rejEtmp;
	
	fprintf('Computing global joint probability...\n');
	tmpdata = permute(icaacttmp, [3 1 2]);
	tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
	[ EEG.stats.icakurt rej] = rejkurt( tmpdata, globthresh, EEG.stats.icakurt, 1); 
end;
rej = rej' | max(rejE, [], 1);
fprintf('%d/%d trials marked for rejection\n', sum(rej), EEG.trials);

if calldisp
	if vistype == 1 % EEGPLOT -------------------------
	    if icacomp == 1 macrorej  = 'EEG.reject.rejkurt';
	        			macrorejE = 'EEG.reject.rejkurtE';
	    else			macrorej  = 'EEG.reject.icarejkurt';
	        			macrorejE = 'EEG.reject.icarejkurtE';
	    end;
		
		colrej = EEG.reject.rejkurtcol;
		eeg_rejmacro; % script macro for generating command and old rejection arrays

	    if icacomp == 1
	        eegplot( EEG.data(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    else
	        eegplot( icaacttmp(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    end;	
    else % REJECTRIALS -------------------------
	  	if icacomp	== 1 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( EEG.data, EEG.stats.kurtE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.kurt, ...
						'threshold', locthresh, 'thresholdg', globthresh, 'normalize', 'off'  );
		else 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( icaacttmp, EEG.stats.icakurtE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.icakurt, ...
						'threshold', locthresh, 'thresholdg', globthresh, 'normalize', 'off' );
		end;		
		nrej = n;
	end;	
else
	% compute rejection locally
	rejtmp = max(rejE(elecrange,:),[],1);
	rej = rejtmp | rej;
	nrej =  sum(rej);
	fprintf('%d trials marked for rejection\n', nrej);
end;
if ~isempty(rej)
	if icacomp	== 1
		EEG.reject.rejkurt = rej;
		EEG.reject.rejkurtE = rejE;
	else
		EEG.reject.icarejkurt = rej;
		EEG.reject.icarejkurtE = rejE;
	end;
end;
nrej = sum(rej);

com = [ com sprintf('%s = pop_rejkurt(%s,%s);', inputname(1), ...
		inputname(1), vararg2str({icacomp,elecrange,locthresh,globthresh,superpose,reject})) ];
if nargin < 3 & nargout == 2
	locthresh = com;
end;

return;
