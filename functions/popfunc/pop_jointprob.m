% pop_jointprob() - rejection of artifact in a dataset using joint 
%                   probability (i.e. probability of activity).
%
% Usage:
%   >> pop_jointprob( INEEG, typerej) % popup
%   >> [OUTEEG, locthresh, globthresh, nrej] = ...
%		= pop_jointprob( INEEG, typerej, electrodes, ...
%                   locthresh, globthresh, superpose, reject, vistype);
%
% Inputs:
%   INEEG      - input dataset
%   typerej    - type of rejection (0 = independent components; 1 = eeg
%              data). Default is 1. For independent components, before
%              thresholding, the activity is normalized for each 
%              component.
%   electrodes - [e1 e2 ...] electrodes (number) to take into 
%              consideration for rejection
%   locthresh  - electrodes probability of activity thresholds in terms
%              of standard-dev.
%   globthresh - global threshold (where all electrode activity are 
%              regrouped).
%   superpose  - 0=do not superpose pre-labelling with previous
%              pre-labelling (stored in the dataset). 1=consider both
%              pre-labelling (using different colors). Default is 0.
%   reject     - 0=do not reject labelled trials (but still store the 
%              labels. 1=reject labelled trials. Default is 0.
%   vistype    - visualization type. 0 is rejepoch() and 1 is eegplot()
%              default is 0.  
%
% Outputs:
%   OUTEEG     - output dataset with updated joint probability array
%   locthresh  - electrodes probability of activity thresholds in terms
%              of standard-dev.
%   globthresh - global threshold (where all electrode activity are 
%              regrouped).
%   nrej       - number of rejected sweeps
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: jointprob(), rejstatepoch(), eegplot(), eeglab(), pop_rejepoch()  

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
% Revision 1.8  2002/08/07 22:33:42  arno
% editing header
%
% Revision 1.7  2002/07/31 00:04:10  arno
% debugging
%
% Revision 1.6  2002/07/30 23:21:55  arno
% debugging
%
% Revision 1.5  2002/07/30 22:27:55  arno
% new type of rejection
%
% Revision 1.4  2002/07/26 17:29:54  arno
% output command debug
%
% Revision 1.3  2002/07/26 17:21:03  arno
% return command if 2 outputs
%
% Revision 1.2  2002/07/26 16:40:54  arno
% switching icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-08-02 add eeglab options -ad

function [EEG, locthresh, globthresh, nrej, com] = pop_jointprob( EEG, icacomp, elecrange, ...
                       		locthresh, globthresh, superpose, reject, vistype, topcommand);
com = '';
if nargin < 1
   help pop_jointprob;
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
               		'Cumulate/compare with previous rejection', ...
         			'Actually reject marked trial(s) (YES or NO)', ...
         			'visualization type (REJECTRIALS|EEGPLOT)' };
	inistr      = { fastif(icacomp, ['1:' int2str(EEG.nbchan)], ['1:' int2str(size(EEG.icaweights,1))])...
					fastif(icacomp, '3', '5'),  ...
					fastif(icacomp, '3', '5'), ...
               		'YES', ...
            		'NO', ...
            		'REJECTRIALS' };

	result       = inputdlg( promptstr, fastif( ~icacomp, 'Reject. improbable comp. -- pop_jointprob()', 'Reject improbable data -- pop_jointprob()'), 1,  inistr);
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
	fprintf('Computing joint probability for channels...\n');
	[ EEG.stats.jpE rejE ] = jointprob( EEG.data, locthresh, EEG.stats.jpE, 1); 

	fprintf('Computing global joint probability...\n');
	tmpdata = permute(EEG.data, [3 1 2]);
	tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
	[ EEG.stats.jp rej ] = jointprob( tmpdata, globthresh, EEG.stats.jp, 1); 
else
	eeg_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
	    tmpdata = EEG.icaact;
	else
        tmpdata = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        tmpdata = reshape( tmpdata, size(tmpdata,1), EEG.pnts, EEG.trials);
    end;
	fprintf('Computing joint probability for components...\n');
    [ EEG.stats.icajpE rejE ] = jointprob( tmpdata, locthresh, EEG.stats.icajpE, 1); 

	fprintf('Computing global joint probability...\n');
	tmpdata2 = permute(tmpdata, [3 1 2]);
	tmpdata2 = reshape(tmpdata2, size(tmpdata2,1), size(tmpdata2,2)*size(tmpdata2,3));
	[ EEG.stats.icajp  rej] = jointprob( tmpdata2, globthresh, EEG.stats.icajp, 1); 
	clear tmpdata2;
end;
rej = rej' | max(rejE, [], 1);
fprintf('%d/%d trials rejected\n', sum(rej), EEG.trials);

if calldisp
	if vistype == 1 % EEGPLOT -------------------------
	    if icacomp == 1 macrorej  = 'EEG.reject.rejjp';
	        			macrorejE = 'EEG.reject.rejjpE';
	    else			macrorej  = 'EEG.reject.icarejjp';
	        			macrorejE = 'EEG.reject.icarejjpE';
	    end;
		colrej = EEG.reject.rejjpcol;
		eeg_rejmacro; % script macro for generating command and old rejection arrays

	    if icacomp == 1
	        eegplot( EEG.data(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    else
	        eegplot( tmpdata(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    end;	
    else % REJECTRIALS -------------------------
	  	if icacomp	== 1 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( EEG.data, EEG.stats.jpE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.jp, ...
						'threshold', locthresh, 'thresholdg', globthresh, 'normalize', 'off'  );
		else 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( tmpdata, EEG.stats.icajpE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.icajp, ...
						'threshold', locthresh, 'thresholdg', globthresh, 'normalize', 'off' );
		end;		
		nrej = n;
	end;	
else
	% compute rejection locally
	rejtmp = max(rejE(elecrange,:),[],1);
	rej = rejtmp | rej;
	nrej =  sum(rej);
	fprintf('%d trials rejected\n', nrej);
end;
if ~isempty(rej)
	if icacomp	== 1
		EEG.reject.rejjp = rej;
		EEG.reject.rejjpE = rejE;
	else
		EEG.reject.icarejjp = rej;
		EEG.reject.icarejjpE = rejE;
	end;
end;
nrej = sum(rej);

com = [ com sprintf('%s = pop_jointprob(%s,%s);', inputname(1), ...
		inputname(1), vararg2str({icacomp,elecrange,locthresh,globthresh,superpose,reject})) ]; 
if nargin < 3 & nargout == 2
	locthresh = com;
end;

return;
