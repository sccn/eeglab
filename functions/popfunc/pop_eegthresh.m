% pop_eegthresh() - reject artifacts by detecting outlier values.  This has 
%                   long been a standard method for selecting data to reject.
%                   Applied either for electrode data or component activations.
% Usage:
%   >> pop_eegthresh( INEEG, typerej); % pop-up interactive window
%   >> [EEG Indexes] = pop_eegthresh( INEEG, typerej, elec_comp, lowthresh, ...
%                upthresh, starttime, endtime, superpose, reject);
%
% Graphic interface:
%   "Electrode|Component number(s)" - [edit box] number(s) of the electrode(s) or 
%                 component(s) to take into consideration. Same as the 'elec_comp'
%                 parameter from the command line.
%   "Lower limit(s)" - [edit box] lower threshold limit(s) (in uV|std. dev.). 
%                 Sets command line parameter 'lowthresh'.
%   "Upper limit(s)" - [edit box] upper threshold limit(s) (in uV|std. dev.). 
%                 Sets command line parameter 'upthresh'.
%   "Start time(s)" - [edit box] starting time limit(s) (in seconds). 
%                 Sets command line parameter 'starttime'.
%   "End time(s)" - [edit box] ending time limit(s) (in seconds). 
%                 Sets command line parameter 'endtime'.
%   "Display with previously marked rejections?" - [edit box] either YES or
%                 NO. This edit box sets command line option 'superpose'.
%   "Reject marked trials?" - [edit box] either YES or NO. This edit box 
%                 sets command line option 'reject'.
%
% Inputs:
%   INEEG      - input EEG dataset
%   typerej    - type of rejection (0 = independent components; 1 = raw
%              data). Default is 1. For independent components, before
%              thresholding the activations are normalized (to have std. dev. 1).
%   elec_comp  - [e1 e2 ...] electrode|component numbers to take 
%              into consideration for rejection
%   lowthresh  - lower threshold limit (in uV|std. dev. For components, the 
%              threshold(s) are in std. dev.). Can be an array if more than one 
%              electrode|component number is given in elec_comp (above). 
%              If fewer values than the number of electrodes|components, the 
%              last value is used for the remaining electrodes|components. 
%   upthresh   - upper threshold limit (in uV|std dev) (see lowthresh above)
%   starttime  - rejection window start time(s) in seconds (see lowthresh above)
%   endtime    - rejection window end time(s) in seconds (see lowthresh)
%   superpose  - [0|1] 0=do not superpose rejection markings on previous
%              rejection marks stored in the dataset: 1=show both current and
%              previously marked rejections using different colors. {Default: 0}.
%   reject     - [1|0] 0=do not actually reject the marked trials (but store the 
%              marks: 1=immediately reject marked trials. {Default: 1}.
% Outputs:
%   Indexes    - index of rejected trials
%     When eegplot() is called, modifications are applied to the current 
%     dataset at the end of the call to eegplot() when the user presses 
%     the 'Reject' button.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegthresh(), eeglab(), eegplot(), pop_rejepoch() 

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

function [EEG, Irej, com] = pop_eegthresh( EEG, icacomp, elecrange, negthresh, posthresh, ...
   						starttime, endtime, superpose, reject, topcommand);

Irej = [];
com = '';
if nargin < 1
   help pop_eegthresh;
   return;
end;  
if nargin < 2
   icacomp = 1;
end;  

if icacomp == 0
	if isempty( EEG.icasphere )
		disp('Error: you must run ICA first'); return;
	end;
end;
if exist('reject') ~= 1
    reject = 1;
end;

if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component (number(s), Ex: 2 4 5):', 'Electrode (number(s), Ex: 2 4 5):'), ...
					fastif(icacomp==0, 'Lower limit(s) (std. dev,  Ex: -3 -2.5 -2):', 'Lower limit(s) (uV, Ex:-20 -10 -15):'), ...
					fastif(icacomp==0, 'Upper limit(s) (std. dev, Ex: 2 2 2.5):', 'Upper limit(s) (uV, Ex: 20 10 15):'), ...
					'Start time(s) (seconds, Ex -0.1 0.3):', ...
					'End time(s) (seconds, Ex 0.2):', ...
                    'Display with previously marked rejections? (YES or NO)', ...
         		    'Reject marked trial(s)? (YES or NO)' };
	inistr      = { fastif(icacomp==1, ['1:' int2str(EEG.nbchan)], '1:5'), ...
					fastif(icacomp==1, '-10', '-20'),  ...
					fastif(icacomp==1, '10', '20'), ...
					num2str(EEG.xmin), ...
					num2str(EEG.xmax), ...
                    'NO', ...
            	    'NO' };

	result       = inputdlg2( promptstr, fastif(icacomp == 0, 'Rejection abnormal comp. values -- pop_eegthresh()', ...
											   'Rejection abnormal elec. values -- pop_eegthresh()'), 1,  inistr, 'pop_eegthresh');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	negthresh    = result{2};
	posthresh    = result{3};
	starttime    = result{4};
	endtime      = result{5};
	switch lower(result{6}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{7}), case 'yes', reject=1; otherwise, reject=0; end;
end;

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	negthresh = eval( [ '[' negthresh ']' ]  );
	posthresh = eval( [ '[' posthresh ']' ]  );
    if isstr(starttime)
        starttime = eval( [ '[' starttime ']' ]  );
    end;
    if isstr(endtime)
        endtime   = eval( [ '[' endtime ']' ]  );
    end;
else
	calldisp = 0;
end;

if any(starttime < EEG.xmin) 
 fprintf('Warning : starttime inferior to minimum time, adjusted\n'); 
	starttime(find(starttime < EEG.xmin)) = EEG.xmin; 
end;
if any(endtime   > EEG.xmax) 
	fprintf('Warning : endtime superior to maximum time, adjusted\n'); 
	endtime(find(endtime > EEG.xmax)) = EEG.xmax;
end;

if icacomp == 1
	[Itmp Irej NS Erejtmp] = eegthresh( EEG.data, EEG.pnts, elecrange, negthresh, posthresh, [EEG.xmin EEG.xmax], starttime, endtime);
    tmpelecIout = zeros(EEG.nbchan, EEG.trials);
    tmpelecIout(elecrange,Irej) = Erejtmp;
else
    icaacttmp = eeg_getdatact(EEG, 'component', elecrange);
	[Itmp Irej NS Erejtmp] = eegthresh( icaacttmp, EEG.pnts, 1:length(elecrange), negthresh, posthresh, [EEG.xmin EEG.xmax], starttime, endtime);
    tmpelecIout = zeros(size(EEG.icaweights,1), EEG.trials);
    tmpelecIout(elecrange,Irej) = Erejtmp;
end;

fprintf('%d channel selected\n', size(elecrange(:), 1));
fprintf('%d/%d trials marked for rejection\n', length(Irej), EEG.trials);
tmprejectelec = zeros( 1, EEG.trials);
tmprejectelec(Irej) = 1;

rej  = tmprejectelec;
rejE = tmpelecIout;
if calldisp
    if icacomp == 1 macrorej  = 'EEG.reject.rejthresh';
        			macrorejE = 'EEG.reject.rejthreshE';
    else			macrorej  = 'EEG.reject.icarejthresh';
        			macrorejE = 'EEG.reject.icarejthreshE';
    end;
	
	colrej = EEG.reject.rejthreshcol;
	eeg_rejmacro; % script macro for generating command and old rejection arrays
	     
    if icacomp == 1
        eegplot( EEG.data(elecrange,:,:), 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
    else
        eegplot( icaacttmp, 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
    end;
else 
    if reject == 1
        EEG = pop_rejepoch(EEG, rej, 0);
    end;
end;
if ~isempty(rej)
    if icacomp	== 1
        EEG.reject.rejthresh  = rej;
        EEG.reject.rejthreshE = rejE;
    else
        EEG.reject.icarejthresh  = rej;
        EEG.reject.icarejthreshE = rejE;
    end;
end;

%com = sprintf('Indexes = pop_eegthresh( %s, %d, [%s], [%s], [%s], [%s], [%s], %d, %d);', ...
%   inputname(1), icacomp, num2str(elecrange),  num2str(negthresh), ...
%   num2str(posthresh), num2str(starttime ) , num2str(endtime), superpose, reject ); 
com = [ com sprintf('%s = pop_eegthresh(%s,%s);', inputname(1), ...
		inputname(1), vararg2str({icacomp,elecrange,negthresh,posthresh,starttime,endtime,superpose,reject})) ]; 
if nargin < 3
	Irej = com;
end;

return;

% reject artifacts in a sequential fashion to save memory (ICA ONLY)
% -------------------------------------------------------
function [Irej, Erej] = thresh( data, elecrange, timerange, negthresh, posthresh, starttime, endtime);
    Irej    = [];
    Erej    = zeros(size(data,1), size(data,2));
    for index = 1:length(elecrange)
       tmpica = data(index,:,:);
       tmpica = reshape(tmpica, 1, size(data,2)*size(data,3));
       
       % perform the rejection
       % ---------------------	
	   tmpica = (tmpica-mean(tmpica,2)*ones(1,size(tmpica,2)))./ (std(tmpica,0,2)*ones(1,size(tmpica,2)));
	   [I1 Itmprej NS Etmprej] = eegthresh( tmpica, size(data,2), 1, negthresh, posthresh, ...
						timerange, starttime, endtime);
 	   Irej = union(Irej, Itmprej);
 	   Erej(elecrange(index),Itmprej) = Etmprej;
	end;

