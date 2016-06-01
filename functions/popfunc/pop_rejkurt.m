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
%   "Electrode|Component" - [edit box] electrodes or components indices to take 
%                 into consideration for rejection. Same as the 'elec_comp'
%                 parameter from the command line.
%   "Single-channel limit|Single-component limit" - [edit box] kurtosis limit
%                  in terms of standard-dev. Same as 'locthresh' command line
%                 parameter.
%   "All-channel limit|All-component limit" - [edit box] kurtosis limit
%                  in terms of standard-dev (all channel regrouped). Same as 
%                 'globthresh' command line parameter.
%   "visualization type" - [popup menu] can be either 'REJECTRIALS'|'EEGPLOT'.
%                 This correspond to the command line input option 'vistype'.
%   "Display with previous rejection(s)" - [checkbox] This checkbox set the
%                 command line input option 'superpose'.          
%   "Reject marked trial(s)" - [checkbox] This checkbox set the command
%                  line input option 'reject'.
% 
% Inputs:
%   INEEG      - Input dataset
%   typerej    - Type of rejection (0 = independent components; 1 = eeg
%              data). Default is 1. For independent components, before
%              thresholding, the activity is normalized for each 
%              component.
%   elec_comp  - [e1 e2 ...] electrodes (number) to take into 
%              consideration for rejection.
%   locthresh  - Activity kurtosis limit in terms of standard-dev.
%   globthresh - Global limit (for all channel). Same unit as above.
%   superpose  - [0] do not superpose pre-labelling with previous
%              pre-labelling (stored in the dataset). [1] consider 
%              both pre-labelling (using different colors). Default is [0].
%   reject     - [0] do not reject labelled trials (but still  
%              store the labels. [1] reject labelled trials. 
%              Default is [1].
%   vistype    - Visualization type. [0] calls rejstatepoch() and [1] calls
%              eegplot() default is [0].When added to the command line
%              call it will not display the plots if the option 'plotflag'
%              is not set.
%   topcommand - [] Deprecated argument , keep to ensure backward compatibility
%   plotflag   - [1,0] [1]Turns plots 'on' from command line, [0] off.
%              (Note for developers: When called from command line 
%              it will make 'calldisp = plotflag') {Default: 0}
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
% 03-08-02 add eeglab options -ad

function [EEG, locthresh, globthresh, nrej, com] = pop_rejkurt( EEG, icacomp, elecrange, ...
                       		locthresh, globthresh, superpose, reject, vistype, topcommand, plotflag); 
nrej = []; com = '';
if nargin < 1
   help pop_rejkurt;
   return;
end;  
if nargin < 2
   icacomp = 1;
end;
if exist('reject') ~= 1
    reject = 1;
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
    promptstr = { [ fastif(icacomp, 'Electrode', 'Component') ' (indices; Ex: 2 6:8 10):' ], ...
                  [ fastif(icacomp, 'Single-channel', 'Single-component') ' limit(s) (std. dev(s): Ex: 2 2 2 2.5):'], ...
                  [ fastif(icacomp, 'All-channel', 'All-component') ' limit(s) (std. dev(s): Ex: 2.1 2 2 2):'], ...
                  'Visualization mode',...
                  'Display previous rejection marks', ...
                  'Reject marked trial(s)'};
    inistr = { fastif(icacomp, ['1:' int2str(EEG.nbchan)], ['1:' int2str(size(EEG.icaweights,1))])...
               fastif(icacomp, '3', '5'),  ...
               fastif(icacomp, '3', '5'), ...
               ''...
               '1', ...
               '0',};
    vismodelist= {'REJECTTRIALS','EEGPLOT'};
    g1  = [1 0.1 0.75];
    g2 = [1 0.26 0.9];
    g3 = [1 0.22 0.85];
    geometry = {g1 g1 g1 g2 [1] g3 g3};
    
    uilist = {...
        { 'Style', 'text', 'string', promptstr{1}} {} { 'Style','edit'      , 'string' ,inistr{1} 'tag' 'cpnum'}...
        { 'Style', 'text', 'string', promptstr{2}} {} { 'Style','edit'      , 'string' ,inistr{2} 'tag' 'singlelimit'}...
        { 'Style', 'text', 'string', promptstr{3}} {} { 'Style','edit'      , 'string' ,inistr{3} 'tag' 'alllimit'}...
        { 'Style', 'text', 'string', promptstr{4}} {} { 'Style','popupmenu' , 'string' , vismodelist 'tag' 'specmethod' }...
        {}...
        { 'Style', 'text', 'string', promptstr{5}} {} { 'Style','checkbox'  ,'string'  ,' ' 'value'  str2double(inistr{5})  'tag' 'rejmarks' }...
        { 'Style', 'text', 'string', promptstr{6}} {} { 'Style','checkbox'  ,'string'  ,' ' 'value'  str2double(inistr{6})  'tag' 'rejtrials'} ...
        };
    figname = fastif(~icacomp, 'Trial rejection using comp. kurtosis -- pop_rejkurt()', 'Trial rejection using data kurtosis -- pop_rejkurt()');
    result = inputgui( geometry,uilist,'pophelp(''pop_rejkurt'');', figname);
   
    size_result  = size( result );
    if size_result(1) == 0, locthresh = []; globthresh = []; return; end;
    elecrange    = result{1};
    locthresh    = result{2};
    globthresh   = result{3};
    switch result{4}, case 1, vistype=0; otherwise, vistype=1; end;
    superpose    = result{5};
    reject       = result{6};
end;

if ~exist('vistype','var'),   vistype   = 0;   end;
if ~exist('reject','var'),    reject    = 0;    end;
if ~exist('superpose','var'), superpose = 1; end;

if isstr(elecrange) % convert arguments if they are in text format 
    calldisp   = 1;
	elecrange  = eval( [ '[' elecrange ']'  ]  );
	locthresh  = eval( [ '[' locthresh ']'  ]  );
	globthresh = eval( [ '[' globthresh ']' ]  );
else
    calldisp = 0;
end;

if exist('plotflag','var') && ismember(plotflag,[1,0])
    calldisp = plotflag;
else
    plotflag = 0;
end

if isempty(elecrange)
	error('No electrode selectionned');
end;	

% compute the joint probability
% -----------------------------
if icacomp == 1
	fprintf('Computing kurtosis for channels...\n');
    tmpdata = eeg_getdatact(EEG);
    if isempty(EEG.stats.kurtE )
		[ EEG.stats.kurtE rejE ] = rejkurt( tmpdata, locthresh, EEG.stats.kurtE, 1); 
	end;
	[ tmp rejEtmp ] = rejkurt( tmpdata(elecrange, :,:), locthresh, EEG.stats.kurtE(elecrange, :), 1); 
    rejE    = zeros(EEG.nbchan, size(rejEtmp,2));
	rejE(elecrange,:) = rejEtmp;
	
	fprintf('Computing all-channel kurtosis...\n');
	tmpdata2 = permute(tmpdata, [3 1 2]);
	tmpdata2 = reshape(tmpdata2, size(tmpdata2,1), size(tmpdata2,2)*size(tmpdata2,3));
	[ EEG.stats.kurt rej ] = rejkurt( tmpdata2, globthresh, EEG.stats.kurt, 1); 
else
	fprintf('Computing joint probability for components...\n');
    % test if ICA was computed
    % ------------------------
    icaacttmp = eeg_getica(EEG);
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
	        eegplot( tmpdata(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    else
	        eegplot( icaacttmp(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    end;	
    else % REJECTRIALS -------------------------
	  	if icacomp	== 1 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( tmpdata, EEG.stats.kurtE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.kurt, ...
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
    if reject
        EEG = pop_rejepoch(EEG, rej, 0);
    end;
end;
nrej = sum(rej);

com = [ com sprintf('%s = pop_rejkurt(%s,%s);', inputname(1), ...
		inputname(1), vararg2str({icacomp,elecrange,locthresh,globthresh,superpose,reject,vistype, [], plotflag})) ];
if nargin < 3 & nargout == 2
	locthresh = com;
end;

return;
