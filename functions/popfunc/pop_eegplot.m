% pop_eegplot() - manual rejection of artifact in a dataset.
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
%              labels. 1=reject labelled trials. Default is 0.
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

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-27-02 added event latency recalculation for continuous data -ad

function [com] = pop_eegplot( EEG, typerej, superpose, reject, topcommand);

com = '';
if nargin < 2
	help pop_eegplot;
	return;
end;	

typerej = ~typerej;
if typerej == 1
	if isempty( EEG.icasphere )
		disp('Error: you must run ICA first'); return;
	end;
end;

superpose=0; 
reject = 1;
if nargin < 3 & EEG.trials > 1

	% which set to save
	% -----------------
    promptstr    = { 'Add to previously labelled rejections (yes/no)', ...
         	     'Reject labelled trials (yes/no)' };
	inistr       = { 'yes', 'no' };
	result       = inputdlg( promptstr, fastif(typerej, 'Manual component rejection -- pop_eegplot()', 'Manual trials rejection -- pop_eegplot()'), 1,  inistr);
	size_result  = size( result );
	if size_result(1) == 0 return; end;
   
   switch lower(result{1}), case 'yes', superpose=1; end;
   switch lower(result{2}), case 'no', reject=0; end;

end;

if EEG.trials > 1
    if typerej == 0 macrorej  = 'EEG.reject.rejmanual';
        			macrorejE = 'EEG.reject.rejmanualE';
    else			macrorej  = 'EEG.reject.icarejmanual';
        			macrorejE = 'EEG.reject.icarejmanualE';
    end;
    elecrange = [1:EEG.nbchan];
	eeg_rejmacro; % script macro for generating command and old rejection arrays

else % case of a single trial (continuous data)
	     %if typerej, 
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
         [  'tmpp = eegplot2event(TMPREJ, -1);' ...
            'if ~isempty(tmpp),' ... 
            '  if isfield(EEG.event, ''latency''), tmpalllatencies = cell2mat( { EEG.event.latency } );' ...
            '  else tmpalllatencies = []; end;' ...
            '  [EEG.data EEG.xmax tmpalllatencies] = eegrej(' fastif(option_keepdataset, 'EEG.data', '''EEG.data'''),  ', tmpp(:,3:4), EEG.xmax-EEG.xmin, tmpalllatencies);' ...
            '  EEG.pnts = size(EEG.data,2);' ...
            '  EEG.xmax = EEG.xmax+EEG.xmin;' ...
            '  if ~isempty(tmpalllatencies)' ...
            '      tmpnanloc = find(~isnan(tmpalllatencies));' ...
            '      EEG.event = EEG.event(tmpnanloc);' ...
            '      fprintf(''Pop_manual: event latencies recomputed and %d event removed (out of %d)\n'', length(tmpalllatencies)-length(EEG.event), length(tmpalllatencies));' ...
            '      tmpalllatencies = tmpalllatencies(tmpnanloc);' ...
            '      for tmpindex = 1:length(EEG.event)' ...
            '          EEG.event = setfield(EEG.event, { tmpindex }, ''latency'', tmpalllatencies(tmpindex));' ...
            '      end;' ...
            '  end;' ...
            '  EEG.icaact = [];' ...
            '  eeg_store; clear tmpalllatencies tmpindex tmpnanloc; eeg_updatemenu;' ...
            'end;' ...
            'clear tmpp;' ];
      %end;
      oldrej = []; oldrejE = [];
end;

if typerej == 0
	eeg_multieegplot( EEG.data, [], [], oldrej, oldrejE, 'title', 'Scroll channel activities -- eegplot()', 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, 'eloc_file', EEG.chanlocs); 
else
	eeg_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
	    tmpdata = EEG.icaact;
	else
        tmpdata = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        tmpdata = reshape( tmpdata, size(tmpdata,1), EEG.pnts, EEG.trials);
    end;
	eeg_multieegplot( tmpdata, [], [], oldrej, oldrejE, 'title', 'Scroll component activities -- eegplot()', 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command); 
end;

com = [ com sprintf('pop_eegplot( %s, %d, %d, %d);', inputname(1), ~typerej, superpose, reject) ]; 
return;
