% pop_rejkurt() - rejection of artifact in a dataset using kurtosis 
%                 of activity (i.e. to detect peaky distribution of
%                 activity).
%
% Usage:
%   >> [OUTEEG, locthresh, globthresh, nrej] = ...
%		= pop_rejkurt( INEEG, typerej, electrodes, ...
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
%   locthresh  - activity kurtosis thresholds in terms of standard-dev.
%   globthresh - global threshold (for all electrodes of a trial).
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
   icacomp = 0;
end;  
icacomp = ~icacomp;
if icacomp == 1
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
	promptstr   = { fastif(icacomp, 'Component (number; ex: 2 4 5):', 'Electrode (number; ex: 2 4 5):'), ...
					'Local thresholds (standard deviation: ex: 3 4 2):', ...
					'Global thresholds (standard deviation: ex: 3 4 2):', ...
               		'Cumulate/compare with current rejection', ...
         			'Actually reject trial (YES or NO for just labelling them', ...
         			'visualization type (REJECTRIALS|EEGPLOT)' };
	inistr      = { ['1:' int2str(EEG.nbchan)], ...
					fastif(icacomp, '5', '3'),  ...
					fastif(icacomp, '5', '3'), ...
               		'YES', ...
            		'NO', ...
            		'REJECTRIALS' };

	result       = inputdlg( promptstr, fastif(icacomp, 'Trial rejection using component kurtosis -- pop_rejkurt()', 'Trial rejection using kurtosis -- pop_rejkurt()'), 1,  inistr);
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
if icacomp == 0
	fprintf('Computing joint probability for channels...\n');
	[ EEG.stats.kurtE rejE ] = rejkurt( EEG.data, locthresh, EEG.stats.kurtE, 1); 

	fprintf('Computing global joint probability...\n');
	tmpdata = permute(EEG.data, [3 1 2]);
	tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
	[ EEG.stats.kurt rej ] = rejkurt( tmpdata, globthresh, EEG.stats.kurt, 1); 
else
	fprintf('Computing joint probability for components...\n');
    % test if ICA was computed
    % ------------------------
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
 	if option_computeica  
    	icaacttmp = EEG.icaact(elecrange, :, :);
	else
        icaacttmp = (EEG.icaweights(elecrange,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        icaacttmp = reshape( icaacttmp, length(elecrange), EEG.pnts, EEG.trials);
    end;
	[ EEG.stats.icakurtE rejE ] = rejkurt( icaacttmp, locthresh, EEG.stats.icakurtE, 1); 

	fprintf('Computing global joint probability...\n');
	tmpdata = permute(icaacttmp, [3 1 2]);
	tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
	[ EEG.stats.icakurt rej] = rejkurt( tmpdata, globthresh, EEG.stats.icakurt, 1); 
end;
rej = rej' | max(rejE, [], 1);
fprintf('%d/%d trials rejected\n', sum(rej), EEG.trials);

if calldisp
	if vistype == 1 % EEGPLOT -------------------------
	    if icacomp == 0 macrorej  = 'EEG.reject.rejkurt';
	        			macrorejE = 'EEG.reject.rejkurtE';
	    else			macrorej  = 'EEG.reject.icarejkurt';
	        			macrorejE = 'EEG.reject.icarejkurtE';
	    end;
		eeg_rejmacro; % script macro for generating command and old rejection arrays

	    if icacomp == 0
	        eeg_multieegplot( EEG.data(elecrange,:,:), rej, rejE, oldrej, oldrejE, 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
	    else
	        eeg_multieegplot( icaacttmp(elecrange,:,:), rej, rejE, oldrej, oldrejE, 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
	    end;	
    else % REJECTRIALS -------------------------
	  	if icacomp	== 0 
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
	fprintf('%d trials rejected\n', nrej);
end;
if ~isempty(rej)
	if icacomp	== 0 
		EEG.reject.rejkurt = rej;
		EEG.reject.rejkurtE = rejE;
	else
		EEG.reject.icarejkurt = rej;
		EEG.reject.icarejkurtE = rejE;
	end;
end;
nrej = sum(rej);

com = [ com sprintf('Indexes = pop_rejkurt( %s, %d, [%s], [%s], [%s], %d, %d);', ...
   inputname(1), ~icacomp, num2str(elecrange),  num2str(locthresh), ...
   num2str(globthresh), superpose, reject ) ]; 

return;
