% eeg_rejmacro() - Internal EEGLAB macro for all pop_ functions that
%                  perform data rejection.   
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

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
% Revision 1.16  2002/08/08 00:23:22  arno
% editing
%
% Revision 1.15  2002/08/08 00:22:36  arno
% adding colmodif option
%
% Revision 1.14  2002/08/07 23:13:16  arno
% updating message
%
% Revision 1.13  2002/08/05 16:49:53  arno
% debugging
%
% Revision 1.12  2002/07/31 17:12:10  arno
% [6~[6~[6~debugging
%
% Revision 1.11  2002/07/31 16:59:04  arno
% debugging
%
% Revision 1.10  2002/07/31 16:38:50  arno
% debugging
%
% Revision 1.9  2002/07/31 16:30:37  arno
% special case for manual rejection
%
% Revision 1.8  2002/07/31 16:16:31  arno
% no change
%
% Revision 1.7  2002/07/31 01:01:57  arno
% debugging for frequencies
%
% Revision 1.6  2002/07/30 23:39:30  arno
% implement rject superposition
%
% Revision 1.5  2002/07/26 16:50:33  arno
% checking icacomp
%
% Revision 1.4  2002/07/08 22:02:42  arno
% adding a warning for data epoch labelling
%
% Revision 1.3  2002/06/25 01:48:19  arno
% changing inputaname to EEG
%
% Revision 1.2  2002/04/26 21:27:36  arno
% updating call to eeg_store
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-08-02 include past rejections in eegplot -ad

% this macro is used by pop_eegplot, pop_eegthresh, pop_const, pop_jointprob, pop_eegthreshfreq, pop_kurtosis

% inputs:
%   rej - input rejection in this workspace (string)
%   rejE - electrode rejection in this workspace (string)
%   elecrange - range of selected electrodes
%   reject - superpose with previous rejections
% optional : nbpnts = number of points per trial

% outputs:
%   rejection arrays ('oldrej', 'oldrejE'), 'command' for eeglot

if ~exist('elecrange')
	help eeg_rejmacro;
	error('Error: eeg_rejmacro can not be called from the command line');
end;	

if ~exist('nbpnts') nbpnts = EEG.pnts; end;
com2 = [ 'if ~isempty(TMPREJ), [tmprej tmprejE] = eegplot2trial(TMPREJ,' ...
		        int2str(nbpnts) ', EEG.trials, [' num2str(colrej) ';' num2str(EEG.reject.rejmanualcol) ...
		 '], []);' ... % include past rejection
         'if ~isempty(tmprejE),' ...
		 '   tmprejE2 = zeros(EEG.nbchan, length(tmprej));' ...
         '   tmprejE2([' int2str(elecrange) '],:) = tmprejE;' ... 
		 'else,' ...
		 '   tmprejE2 = [];' ...
		 'end;' ...
         macrorej '= tmprej;' macrorejE '= tmprejE2;' ];
if reject
    com2 = [com2 sprintf(['[%s tmpcom] = pop_rejepoch(%s, tmprej, 1);' ...
		   'if ~isempty(tmpcom),'...
		   '   [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET);' ...
		   '    disp(''warning: manual modification were not saved in the history'');' ...
	       'end; eeglab(''redraw''); end;'], 'EEG', 'EEG'); ] ;
else
	com2 = [com2 ...
		    'disp(''warning: manual modification were not saved in the history'');' ...
			' warndlg2(strvcat(''Epochs labeled for rejection have been noted'',' ...
			'''To actually reject these epochs, use'', ''Tools/Reject data epochs/Reject labeled epochs''), ''Warning'');' ...
			'[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab(''redraw''); end;' ];
end; 
if ~exist('topcommand')
	topcommand = [];
end;

% the first part is used to convert the eegplot output
command = [  com2 topcommand 'clear tmpcom tmprej tmprejE tmprejE2 TMPREJ;' ];

if all(colrej == EEG.reject.rejmanualcol)
	oldrej = [];  % for manual rejection, old rejection are
	oldrejE = []; % the current rejection
else
	oldrej  = eval(macrorej);
	oldrejE = eval(macrorejE);
end;

% mix all type of rejections
% --------------------------
switch superpose
 case 0, rejeegplot = trial2eegplot(  rej, rejE, nbpnts, colrej);
 case 1, rejeegplottmp = trial2eegplot(  oldrej, oldrejE, nbpnts, min(colrej+0.15, [1 1 1]));
         if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplottmp ]; 
		 else rejeegplot = []; end;
         rejeegplottmp = trial2eegplot(  rej, rejE, nbpnts, colrej);
         if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplot; rejeegplottmp ]; end;
 case 2, 
  rejeegplot = [];
  for index = 1:length(EEG.reject.disprej)
	  if ~isempty(EEG.reject.disprej{index})
		  eval([ 'colortmp = EEG.reject.rej' EEG.reject.disprej{index} 'col;']); 
		  if any(colortmp ~= colrej) % test if current rejection (if color different)
			  if icacomp == 0 % ica
				  currentname = [ 'EEG.reject.icarej' EEG.reject.disprej{index} ];
			  else
				  currentname = [ 'EEG.reject.rej' EEG.reject.disprej{index} ];
			  end;
			  currentcolor =  [ 'EEG.reject.rej' EEG.reject.disprej{index} 'col' ];
			  if strcmp(EEG.reject.disprej{index}, 'manual')
				  currentcolor = [ 'min(' currentcolor '+0.15, [1 1 1])' ];
			  end; % using this test, manual rejections won't be added to current rej
			  eval( [ 'rejeegplottmp = trial2eegplot( ' currentname ',' currentname ...
					  'E, nbpnts,' currentcolor ');' ]);
			  if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplot; rejeegplottmp ]; end;
		  end;
	  end;
  end;
  rejeegplottmp = trial2eegplot(  rej, rejE, nbpnts, colrej);
  if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplot; rejeegplottmp ]; end;
end;
if ~isempty(rejeegplot)
	rejeegplot = rejeegplot(:,[1:5,elecrange]);
else
	rejeegplot = [];
end;
eegplotoptions = { 'winlength', 5, 'position', [100 300 800 500], 'winrej', ...
				   rejeegplot, 'xgrid', 'off', 'wincolor', EEG.reject.rejmanualcol, ...
				   'colmodif', { { EEG.reject.rejmanualcol colrej } } };

if ~isempty(EEG.chanlocs) & icacomp == 1
	eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs };
end;
if ~reject
	eegplotoptions = { eegplotoptions{:}  'butlabel', 'UPDATE MARKS' };
end;
