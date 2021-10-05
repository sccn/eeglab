% pop_rejtrend() - Measure linear trends in EEG data; reject data epochs 
%                  containing strong trends.
% Usage:
%   >> pop_rejtrend( INEEG, typerej); % pop up an interactive window
%   >> OUTEEG = pop_rejtrend( INEEG, typerej, elec_comp, ...
%                winsize, maxslope, minR, superpose, reject,calldisp);
%
% Pop-up window interface:
%   "Electrode|Component indices(s)" - [edit box] electrode or component indices 
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
%   "Display with previous rejection(s)" - [checkbox] This checkbox set the
%                 command line input option 'superpose'.          
%   "Reject marked trial(s)" - [checkbox] This checkbox set the command
%                  line input option 'reject'.
% Command line inputs:
%   INEEG      - input EEG dataset
%   typerej    - [1|0] data to reject on: 0 = component activations; 
%                1 = electrode data. {Default: 1}.
%   elec_comp  - [e1 e2 ...] electrode|component number(s) to take into 
%                consideration during rejection
%   winsize    - (integer) number of consecutive points
%                to use in detecting linear trends
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
% 03-07-02 added srate argument to eegplot call -ad
% 03-07-02 add the eeglab options -ad

function [EEG, com] = pop_rejtrend( EEG, icacomp, elecrange, winsize, ...
				minslope, minstd, superpose, reject, calldisp);

com = '';
if nargin < 1
   help pop_rejtrend
   return;
end
if nargin < 2
   icacomp = 1;
end
if icacomp == 0
	if isempty( EEG.icasphere )
	    ButtonName=questdlg( 'Do you want to run ICA now ?', ...
                         'Confirmation', 'NO', 'YES', 'YES');
    	switch ButtonName
        	case 'NO', disp('Operation cancelled'); return;   
        	case 'YES', [ EEG, com ] = pop_runica(EEG);
    	end % switch
	end
end
if exist('reject') ~= 1
    reject = 1;
end
if nargin < 3
    
    % which set to save
    % -----------------
    promptstr = { [fastif(icacomp, 'Electrode', 'Component') ' (indices; Ex: 2 6:8 10):' ], ...
                  'Slope window width (in points)', ...
                  [fastif(icacomp,'Maximum slope to allow (std. dev./epoch)','Maximum slope to allow (uV/epoch)')], ...
                  'R-square limit to allow ([0:1], Ex: 0.8)', ...
                  'Display previous rejection marks', ...
                  'Reject marked trial(s)' };
    
    inistr = { ['1:' int2str(EEG.nbchan)], ...
               int2str(EEG.pnts),  ...
               '0.5', ...
               '0.3', ...
               '0',   ...
               '0' };
    g1 = [1 0.1 0.75];
    g2 = [1 0.22 0.85];
    geometry = {g1 g1 g1 g1 1 g2 g2};
    
    uilist = {...
              { 'Style', 'text', 'string', promptstr{1}} {} { 'Style','edit'      ,'string' ,inistr{1}  'tag' 'cpnum'}...
              { 'Style', 'text', 'string', promptstr{2}} {} { 'Style','edit'      ,'string' ,inistr{2}  'tag' 'win' }...
              { 'Style', 'text', 'string', promptstr{3}} {} { 'Style','edit'      ,'string' ,inistr{3}  'tag' 'maxslope'}...
              { 'Style', 'text', 'string', promptstr{4}} {} { 'Style','edit'      ,'string' ,inistr{4}  'tag' 'rlim'}...
              {}...
              { 'Style', 'text', 'string', promptstr{5}} {} { 'Style','checkbox'  ,'string'  ,' ' 'value' str2double(inistr{6}) 'tag','rejmarks' }...
              { 'Style', 'text', 'string', promptstr{6}} {} { 'Style','checkbox'  ,'string'  ,' ' 'value' str2double(inistr{6}) 'tag' 'rejtrials'} ...
              };
    
    figname = fastif(~icacomp, 'Trend rejection in component(s) -- pop_rejtrend()','Data trend rejection -- pop_rejtrend()');
    result = inputgui( geometry,uilist,'pophelp(''pop_rejtrend'');', figname);
    size_result  = size( result );
    if size_result(1) == 0 return; end
    elecrange    = result{1};
    winsize      = result{2};
    minslope     = result{3};
    minstd       = result{4};
    superpose    = result{5};
    reject       = result{6};
    calldisp     = 1;
end

calldisp = 0;
if ~exist('superpose','var'), superpose = 0; end
if ~exist('reject','var'),    reject    = 0; end
if ~exist('calldisp','var'),  calldisp  = 1; end

if nargin < 8
    calldisp = 1;
end

if ischar(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	winsize   = eval( [ '[' winsize ']' ]  );
	minslope  = eval( [ '[' minslope ']' ]  );
	minstd    = eval( [ '[' minstd ']' ]  );
end

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
    end
	colrej = EEG.reject.rejconstcol;
	eeg_rejmacro; % script macro for generating command and old rejection arrays

	if icacomp == 1
		eegplot( EEG.data(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	else
		eegplot( icaacttmp, 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
    end
else
    if ~isempty(rej)
        if icacomp	== 1
            EEG.reject.rejconst = rej;
            EEG.reject.rejconstE = rejE;
        else
            EEG.reject.icarejconst = rej;
            EEG.reject.icarejconstE = rejE;
        end
        if reject
            EEG = pop_rejepoch(EEG, rej, 0);
        end
    end
end

com = [ com sprintf('EEG = pop_rejtrend(EEG,%s);', ...
		vararg2str({icacomp,elecrange,winsize,minslope,minstd,superpose, ~calldisp & reject })) ]; 

return;
