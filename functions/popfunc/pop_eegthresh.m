% pop_eegthresh() - rejection of artifact by detecting abnormal values
%               (i.e. standard) method.
%
% Usage:
%   >> pop_eegthresh( INEEG, typerej); % pops-up
%   >> [EEG Indexes] = pop_eegthresh( INEEG, typerej, elec_comp, negthresh, ...
%                posthresh, starttime, endtime, superpose, reject);
%
% Inputs:
%   INEEG      - input dataset
%   typerej    - type of rejection (0 = independent components; 1 = eeg
%              data). Default is 1. For independent components, before
%              thresholding, the activity is normalized.
%   elec_comp  - [e1 e2 ...] electrodes (number) or components to take 
%              into consideration for rejection
%   negthresh  - negative threshold limit in mV (can be an array if 
%              several electrodes; if less numbe  of values than number 
%              of electrodes the last value is used for the remaining 
%              electrodes). For independent component, this threshold is
%              expressed in term of standard deviation. 
%   posthresh  - positive threshold limit in mV (same syntax as negthresh)
%   starttime  - starting time limit in second (same syntax as negthresh)
%   endtime    - ending time limit in second (same syntax as negthresh)
%   superpose  - 0=do not superpose pre-labelling with previous
%              pre-labelling (stored in the dataset). 1=consider both
%              pre-labelling (using different colors). Default is 0.
%   reject     - 0=do not reject labelled trials (but still store the 
%              labels. 1=reject labelled trials. Default is 0.
%
% Outputs:
%   Indexes    - index of rejected sweeps
%     when eegplot is called, modifications are applied to the current 
%     dataset at the end of the call of eegplot (when the user press the 
%     button 'reject').
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegthresh(), eeglab(), eegplot(), pop_rejepoch() 

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
% Revision 1.8  2002/08/07 22:40:09  arno
% same
%
% Revision 1.7  2002/08/07 22:26:58  arno
% editing header
%
% Revision 1.6  2002/07/30 23:59:35  arno
% debugging for ICA
%
% Revision 1.5  2002/07/30 23:26:09  arno
% new rejection type
%
% Revision 1.4  2002/07/26 17:56:04  arno
% still debugging
%
% Revision 1.3  2002/07/26 17:52:24  arno
% debugging
%
% Revision 1.2  2002/07/26 16:48:32  arno
% switching icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

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

if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component (number; ex: 2 4 5):', 'Electrode (number; ex: 2 4 5):'), ...
					fastif(icacomp==0, 'Lower limit(s) (std: ex:-3 -4 -2):', 'Lower limit(s) (uV; ex:-20 -10 -15):'), ...
					fastif(icacomp==0, 'Upper limit(s) (std: ex:3 4 2):', 'Upper limit (uV; ex:20 10 15):'), ...
					'Start time(s) (s;ex -0.1 0.3):', ...
					'End time(s) (s;ex 0.2):', ...
                    'Display with previous rejection', ...
         		    'Actually reject marked trial(s) (YES or NO)' };
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
	starttime = eval( [ '[' starttime ']' ]  );
	endtime   = eval( [ '[' endtime ']' ]  );
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
    % test if ICA was computed
    % ------------------------
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
 	if option_computeica  
    	icaacttmp = EEG.icaact(elecrange, :, :);
	else
        icaacttmp = (EEG.icaweights(elecrange,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        icaacttmp = reshape( icaacttmp, length(elecrange), EEG.pnts, EEG.trials);
    end;
	[Itmp Irej NS Erejtmp] = eegthresh( icaacttmp, EEG.pnts, elecrange, negthresh, posthresh, [EEG.xmin EEG.xmax], starttime, endtime);
    tmpelecIout = zeros(size(EEG.icaweights,1), EEG.trials);
    tmpelecIout(elecrange,Irej) = Erejtmp;
end;

fprintf('%d channel selected\n', size(elecrange(:), 1));
fprintf('%d/%d trials rejected\n', length(Irej), EEG.trials);
tmprejectelec = zeros( 1, EEG.trials);
tmprejectelec(Irej) = 1;

if calldisp
    if icacomp == 1 macrorej  = 'EEG.reject.rejthresh';
        			macrorejE = 'EEG.reject.rejthreshE';
    else			macrorej  = 'EEG.reject.icarejthresh';
        			macrorejE = 'EEG.reject.icarejthreshE';
    end;
	
	colrej = EEG.reject.rejthreshcol;
	rej  = tmprejectelec;
	rejE = tmpelecIout;
	eeg_rejmacro; % script macro for generating command and old rejection arrays
	     
    if icacomp == 1
        eegplot( EEG.data(elecrange,:,:), 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
    else
        eegplot( icaacttmp, 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
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

