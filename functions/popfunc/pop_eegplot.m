% pop_eegplot() - reject by visual inspection of artifact in a dataset.
%
% Usage:
%   >> pop_eegplot( INEEG, typerej, superpose, reject );
%
% Inputs:
%   INEEG      - input dataset
%   typerej    - type of rejection (0 = independent components; 1 = eeg
%              data). Default is 1.
%   superpose  - 0 = do not superpose pre-labelling with previous
%              pre-labelling (stored in the dataset). 1=consider both
%              pre-labelling (using different colors). Default is 0.
%   reject     - 0 = do not reject labelled trials (but still store the 
%              them. 1=reject labelled trials). Default is 0.
%
% Outputs:
%   Modifications are applied to the current dataset at the end of the
%   call of eegplot (when the user press the button 'reject').
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), eegplot(), pop_rejepoch()

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
% Revision 1.14  2002/08/12 16:26:20  arno
% inputdlg2
%
% Revision 1.13  2002/08/08 01:33:56  arno
% adding history for continuous rejection
%
% Revision 1.12  2002/08/07 22:41:01  arno
% editing text
%
% Revision 1.11  2002/07/31 17:16:15  arno
% debugging
%
% Revision 1.10  2002/07/31 16:32:09  arno
% debugging
%
% Revision 1.9  2002/07/30 22:21:11  arno
% debugging
%
% Revision 1.8  2002/07/30 22:19:24  arno
% removing multieegplot.m
%
% Revision 1.7  2002/07/27 01:35:55  arno
% editing header
%
% Revision 1.6  2002/07/26 16:53:22  arno
% switching icacomp
%
% Revision 1.5  2002/07/08 22:03:02  arno
% nothing
%
% Revision 1.4  2002/06/25 22:23:56  scott
% *** empty log message ***
%
% Revision 1.3  2002/04/26 21:28:17  arno
% eeg_updatemenu
%
% Revision 1.2  2002/04/26 21:25:03  arno
% update call to eeg_store
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-27-02 added event latency recalculation for continuous data -ad

function [com] = pop_eegplot( EEG, icacomp, superpose, reject, topcommand);

com = '';
if nargin < 1
	help pop_eegplot;
	return;
end;	
if nargin < 2
	icacomp = 1;
end;	
if nargin < 3
	superpose = 0;
end;
if nargin < 4
	reject = 1;
end;
if icacomp == 0
	if isempty( EEG.icasphere )
		disp('Error: you must run ICA first'); return;
	end;
end;

if nargin < 3 & EEG.trials > 1

	% which set to save
	% -----------------
    promptstr    = { 'Add to previously labelled rejections (yes/no)', ...
         	         'Reject labelled trials (yes/no)', ...
						 };
	inistr       = { 'yes', 'no' };
	result       = inputdlg2( promptstr, fastif(icacomp==0, 'Manual component rejection -- pop_eegplot()', ...
								'Reject epochs by visual inspection -- pop_eegplot()'), 1,  inistr, 'pop_eegplot');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
   
   switch lower(result{1}), case 'yes', superpose=1; end;
   switch lower(result{2}), case 'no', reject=0; end;

end;

if EEG.trials > 1
    if icacomp == 1 macrorej  = 'EEG.reject.rejmanual';
        			macrorejE = 'EEG.reject.rejmanualE';
    else			macrorej  = 'EEG.reject.icarejmanual';
        			macrorejE = 'EEG.reject.icarejmanualE';
    end;
    elecrange = [1:EEG.nbchan];
	colrej = EEG.reject.rejmanualcol;
	rej  = eval(macrorej);
	rejE = eval(macrorejE);
	
	eeg_rejmacro; % script macro for generating command and old rejection arrays

else % case of a single trial (continuous data)
	     %if icacomp, 
    	 %    	command = ['if isempty(EEG.event) EEG.event = [eegplot2event(TMPREJ, -1)];' ...
         %         'else EEG.event = [EEG.event(find(EEG.event(:,1) ~= -1),:); eegplot2event(TMPREJ, -1, [], [0.8 1 0.8])];' ...
         %         'end;']; 
      	 %else, command = ['if isempty(EEG.event) EEG.event = [eegplot2event(TMPREJ, -1)];' ...
         %         'else EEG.event = [EEG.event(find(EEG.event(:,1) ~= -2),:); eegplot2event(TMPREJ, -1, [], [0.8 1 0.8])];' ...
         %         'end;']; 
      	 %end;
         %if reject
         %   command = ...
         %   [  command ...
         %      '[EEG.data EEG.xmax] = eegrej(EEG.data, EEG.event(find(EEG.event(:,1) < 0),3:end), EEG.xmax-EEG.xmin);' ...
         %      'EEG.xmax = EEG.xmax+EEG.xmin;' ...
         %   	'EEG.event = EEG.event(find(EEG.event(:,1) >= 0),:);' ...
         %      'EEG.icaact = [];' ...
         %      'EEG = eeg_checkset(EEG);' ];
         eeg_options; % changed from eeglaboptions 3/30/02 -sm
         command = ...
         [  '[EEGTMP LASTCOM] = eeg_eegrej(EEG,eegplot2event(TMPREJ, -1));' ...
            'if ~isempty(LASTCOM),' ... 
            '  h(LASTCOM);' ...
		    '  [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
            '  if ~isempty(LASTCOM),' ... 
            '     h(LASTCOM);' ...
			'     eeglab(''redraw'');' ...
            '  end;' ...
			'end;' ...
            'clear EEGTMP;' ];
      %end;
	  eegplotoptions = { 'winlength', 5, 'position', [100 300 800 500] };
	  if ~isempty(EEG.chanlocs)
		  eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs };
	  end;
end;

if icacomp == 1
	eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Scroll channel activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	%eeg_multieegplot( EEG.data, [], [], oldrej, oldrejE, 'title', 'Scroll channel activities -- eegplot()', 'srate', ...
	%	      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
else
	eeg_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
	    tmpdata = EEG.icaact;
	else
        tmpdata = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        tmpdata = reshape( tmpdata, size(tmpdata,1), EEG.pnts, EEG.trials);
    end;
	eegplot( tmpdata, 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	%eeg_multieegplot( tmpdata, [], [], oldrej, oldrejE, 'title', 'Scroll component activities -- eegplot()', 'srate', ...
	%	      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
end;

com = [ com sprintf('pop_eegplot( %s, %d, %d, %d);', inputname(1), icacomp, superpose, reject) ]; 
return;
